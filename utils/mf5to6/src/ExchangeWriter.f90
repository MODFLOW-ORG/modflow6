module ExchangeWriterModule
  
  use ConstantsModule, only: MAXCHARLEN
  use ExchangeModule, only: ExchangeType
  use FileWriterModule, only: FileWriterType
  use InputOutputModule, only: GetUnit
  use LGRMODULE, ONLY: ISCHILD,NGRDS,NPLBEG,NPRBEG,NPCBEG,NPLEND, &
                       NPREND,NPCEND,NCPP,NPL,IBOTFLG,ISHFLG,IBFLG, &
                       IUPBHSV,IUCBHSV,IUPBFSV,IUCBFSV,MXLGRITER, &
                       IOUTLGR,NBNODES,NPBNODES,IBMAXH,IBMAXF, &
                       NCMAXH,NCMAXF,RELAXH,RELAXF,HCLOSELGR, &
                       FCLOSELGR,HDIFFM,FDIFFM,PRATIN,CRATIN, &
                       PRATOUT,CRATOUT,IBPFLG,IEDG,JEDG,NCPPL, &
                       NODEH,NODEF,NCON,KPLC,IPLC,JPLC,IFACEGN, &
                       ICBOUND,GNHEAD,DHGN,GNFLUX,GNFLUXR, &
                       GNFLUXOLD,HOLDC,GNCOND,VCB,HK, &
                       VK,LGRDAT, SGWF2LGR2PNT
  use ListModule, only: ListType
  use ModelConverterModule, only: ModelConverterType
  use ModelModule, only: ModelType
  use MvrPackageWriterModule, only: MvrPackageWriterType
  use SimPHMFModule, only: store_error, ustop
  
  type, extends(FileWriterType) :: ExchangeWriterType
    type(ListType) :: Exchanges
    type(MvrPackageWriterType), pointer :: MvrWriter => null()
  contains
    procedure, public :: DefineExchanges
    procedure, public :: InitializeExchangeWriter
    procedure, public :: WriteExchangeFiles
  end type ExchangeWriterType
  
contains

  subroutine InitializeExchangeWriter(this, MvrPackWriter)
    implicit none
    ! dummy
    class(ExchangeWriterType) :: this
    type(MvrPackageWriterType), pointer :: MvrPackWriter
    !
    this%MvrWriter => MvrPackWriter
    !
    return
  end subroutine InitializeExchangeWriter

  subroutine DefineExchanges(this, modelConverters)
    implicit none
    ! dummy
    class(ExchangeWriterType) :: this
    type(ListType) :: modelConverters
    ! local
    integer :: i, ngrids, igridchild, igridparent
    character(len=MAXCHARLEN) :: exgfilename, gncfilename
    type(ModelConverterType), pointer :: childConverter
    type(ModelConverterType), pointer :: parentConverter
    type(ExchangeType), pointer :: exchange
    type(ModelType), pointer :: parent => null(), child => null()
    class(*), pointer :: obj => null()
    ! formats
    10 format()
    !
    ngrids = modelConverters%Count()
    !
    ! Iterate over children; a parent/child
    ! exchange is defined for each child.
    do i=2,ngrids
      ! Get child
      obj => modelConverters%GetItem(i)
      select type (obj)
      type is (ModelConverterType)
        childConverter => obj
      end select
      child => childConverter%model
      igridchild = child%IGrid
      !
      ! Get parent of child
      parentConverter => childConverter%ParentConverter
      parent => parentConverter%model
      igridparent = parent%IGrid
      !
      ! Define parent/child exchange
      allocate(exchange)
      ! Define exchange file name and initialize file
      exgfilename = trim(parent%ModelName) // '_' //  &
                    trim(child%ModelName) // '.exg6'
      exchange%FileWriter%Active = .true.
      call exchange%FileWriter%InitializeFile(exgfilename, 'EXG')
      gncfilename = trim(parent%ModelName) // '_' //  &
                    trim(child%ModelName) // '.gnc6'
      exchange%GncFileWriter%Active = .true.
      call exchange%GncFileWriter%InitializeFile(gncfilename, 'GNC')
      !
      call exchange%InitializeExchange(parentConverter%model, &
                                       childConverter%model)
      ! Add exchange to list of exchanges.
      obj => exchange
      call this%Exchanges%Add(obj)
    enddo
    !
    return
  end subroutine DefineExchanges
  
  subroutine WriteExchangeFiles(this)
    implicit none
    ! dummy
    class(ExchangeWriterType) :: this
    ! local
    integer :: i, nexg
    class(*), pointer :: obj => null()
    type(ExchangeType), pointer :: exg => null()
    !
    nexg = this%Exchanges%Count()
    do i=1,nexg
      obj => this%Exchanges%GetItem(i)
      select type (obj)
      type is (ExchangeType)
        exg => obj
        call exg%WriteExchangeFile(this%MvrWriter%fileobj%FName)
        call exg%WriteGhostNodeFile()
      end select
    enddo
    !
    return
  end subroutine WriteExchangeFiles
  
end module ExchangeWriterModule
