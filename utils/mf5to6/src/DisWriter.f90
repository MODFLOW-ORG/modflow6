module DisWriterModule

  use ConstantsModule, only: MAXCHARLEN
  use ConstantsPHMFModule, only: HUGEDBL
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: NCOL, NROW, NLAY, NBOTM, IXSEC, LBOTM, LAYCBD, &
                    DELR, DELC, BOTM, constantdelr, constantdelc, &
                    IBOUND, NCNFBD
  use GlobalVariablesModule, only: echo
  use GWFBCFMODULE, only: LAYCON
  use SimPHMFModule, only: store_error, store_note, store_warning, ustop
  use UtilitiesModule, only: ConstantInt2D, ConstantInt3D, &
                             ConstantReal1D, ConstantReal2D, &
                             ConstantReal3D, Write2Dint, Write3Dint, &
                             Write1Drel, Write2Drel, Write3Drel

  private
  public :: DisWriterType

  type, extends(FileWriterType) :: DisWriterType
    logical, public :: ConstantDelr
    logical, public :: ConstantDelc
    logical, public :: ConstantTop
    integer, public :: Nlaytot
    integer,                   pointer, public  :: Nlaynew => null()
    integer, dimension(:),     pointer, public  :: Layptr => null()
    character(len=12),                  private :: LengthUnit = 'UNDEFINED'
    double precision, dimension(:,:), pointer, private :: TopNew => null()
    double precision :: TopMax = -HUGEDBL
    double precision :: BotMin = HUGEDBL
  contains
    ! Public procedures
    procedure, public :: InitializeFile => initialize_dis
    procedure, public :: WriteDisFile
    procedure, public :: AssignLengthUnit
    ! Private procedures
    procedure, private :: DefineTopNew
    procedure, private :: WriteOptions
    procedure, private :: WriteDimensions
    procedure, private :: WriteDisData
    procedure, private :: WriteIdomain
  end type DisWriterType

contains

  subroutine initialize_dis(this, fname, ftype, pkgname)
    implicit none
    ! dummy
    class(DisWriterType), intent(inout) :: this
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: ftype
    character(len=*), intent(in), optional :: pkgname
    !
    ! Invoke superclass initializer
    if (present(pkgname)) then
      call this%FileWriterType%InitializeFile(fname, ftype, pkgname)
    else
      call this%FileWriterType%InitializeFile(fname, ftype)
    endif
    !
    ! DIS is always active
    this%Active = .true.
    !
    return
  end subroutine initialize_dis

  subroutine WriteDisFile(this, nc, nr, nl, idomain)
    implicit none
    ! dummy
    class(DisWriterType), intent(inout) :: this
    integer, intent(in) :: nc, nr, nl
    integer, pointer, dimension(:,:,:), intent(inout) :: idomain
    ! local
    !
    call this%DefineTopNew()
    call this%WriteOptions()
    call this%WriteDimensions()
    call this%WriteDisData(nc, nr, nl, idomain)
    call this%CloseFile()
    !
    return
  end subroutine WriteDisFile

  subroutine WriteOptions(this)
    implicit none
    ! dummy
    class(DisWriterType), intent(inout) :: this
    ! local
    integer :: iu
    ! formats
    5  format()
    10 format('BEGIN Options')
    20 format(2x,'LENGTH_UNITS',2x,a)
    30 format('END Options')
    !
    iu = this%fileobj%IUnit
    !
    write(iu,5)
    write(iu,10)
    write(iu,20) trim(this%LengthUnit)
    write(iu,30)
    !
    return
  end subroutine WriteOptions

  subroutine WriteDimensions(this)
    implicit none
    ! dummy
    class(DisWriterType), intent(inout) :: this
    ! local
    integer :: iu
    ! formats
    5  format()
    10 format('BEGIN Dimensions')
    20 format('END Dimensions')
    30 format(2x,a,2x,i0)
    !
    iu = this%fileobj%IUnit
    this%Nlaytot = NLAY + NCNFBD
    !
    ! Nlay, Nrow, and Ncol are read in  > SGWF2BAS7ARDIS
    write(iu,5)
    write(iu,10)
    write(iu,30)'NLAY', this%Nlaytot
    write(iu,30)'NROW', nrow
    write(iu,30)'NCOL', ncol
    write(iu,20)
    write(iu,5)
    !
    return
  end subroutine WriteDimensions

  subroutine WriteDisData(this, nc, nr, nl, idomain)
    implicit none
    ! dummy
    class(DisWriterType), intent(inout) :: this
    integer, intent(in) :: nc, nr, nl
    integer, pointer, dimension(:,:,:), intent(in) :: idomain
    ! local
    integer :: iprn, iu, k
    logical :: constant
    double precision :: rval
    character(len=100) :: msg
    ! formats
    10 format(a)
    20 format(2x,a)
    30 format(a,1x,i0)
    !
    iu = this%fileobj%IUnit
    write(iu,10)'BEGIN GRIDDATA'
    !
    ! Write DELR array
    if (echo) then
      iprn = 3
    else
      iprn = -3
    endif
    call ConstantReal1D(ncol, DELR, constant, rval)
    call Write1Drel(iu, ncol, DELR, constant, rval, &
                    'DELR', .true., iprn)
    !
    ! Write DELC array
    call ConstantReal1D(NROW, DELC, constant, rval)
    call Write1Drel(iu, NROW, DELC, constant, rval, &
                    'DELC', .true., iprn)
    !
    ! Write TOP array
    write(iu,20)'TOP LAYERED'
    write(msg,10)'TOP array'
    call ConstantReal2D(NCOL, NROW, this%TopNew, constant, rval)
    call Write2Drel(iu, nrow, ncol, this%TopNew, constant, &
                    rval, msg, .false., iprn)
    !
    ! Write BOTM arrays
    write(iu,20)'BOTM LAYERED'
    do k=1,this%Nlaytot
      write(msg,30)'BOTM layer',k
      call ConstantReal2D(NCOL, NROW, BOTM(:,:,k), constant, rval)
      call Write2Drel(iu, nrow, ncol, BOTM(:,:,k), constant, rval, &
                      msg, .false., iprn)
    enddo
    !
    if (.not. associated(idomain)) then
      call store_error('idomain not associated in DisWriter%WriteDisData')
      call ustop()
    endif
    ! Write Idomain array
    call this%WriteIdomain(nc, nr, nl, idomain)
    !
    write(iu,10)'END GRIDDATA'
    return
  end subroutine WriteDisData

  subroutine WriteIdomain(this, nc, nr, nl, idomain)
    implicit none
    ! dummy
    class(DisWriterType), intent(inout) :: this
    integer, intent(in) :: nc, nr, nl
    integer, dimension(nc,nr,nl), intent(in) :: idomain
    ! local
    integer :: iu, ival, k
    logical :: constant, needToWrite
    character(len=300) :: msg
    ! formats
    20 format(2x,a)
    30 format(a,1x,i0)
    35 format(4x,a,2x,i0)   
    40 format('Idomain for MODFLOW 6 model layer ',i0, &
        ', which replaces a quasi-3D unit, has been assigned equal' &
        ' to Ibound for layer ',i0,' of original model.')
    !
    ! Write IDOMAIN arrays
    call ConstantInt3D(nc, nr, nl, Idomain, constant, ival)
    needToWrite = .true.
    ! IDomain by default = 1, so if ibound==1, no need to write
    if (constant) then
      if (ival==1) needToWrite = .false.
    endif
    !
    if (needToWrite) then
      iu = this%fileobj%IUnit
      write(iu,20)'IDOMAIN LAYERED'
      if (constant) then
        write(iu,35)'CONSTANT',ival
      else
        do k=1,this%Nlaynew
          write(msg,30)'Idomain layer',k
          call Write2Dint(iu, nr, nc, Idomain(:,:,k), constant, &
                          ival, msg, .false., 5)
        enddo
      endif
    endif
    !
    return
  end subroutine WriteIdomain

  subroutine AssignLengthUnit(this, LengthUnit)
    implicit none
    class(DisWriterType), intent(inout) :: this
    character(len=*), intent(in) :: LengthUnit
    !
    this%LengthUnit = LengthUnit
    return
  end subroutine AssignLengthUnit

  subroutine DefineTopNew(this)
    implicit none
    ! dummy
    class(DisWriterType) :: this
    ! local
    integer :: i, j, kbot
    character(len=MAXCHARLEN) :: msg
    ! format
    10 format(g10.3)
    !
    ! Allocate and populate TopNew, and assign this%TopMax
    allocate(this%TopNew(NCOL, NROW))
    do i=1,nrow
      do j=1,ncol
        this%TopNew(j,i) = botm(j,i,0)
        if (botm(j,i,0) > this%TopMax) this%TopMax = botm(j,i,0)
      enddo
    enddo
    !
    ! Find lowest cell bottom
    kbot = LBOTM(NLAY)
    do i=1,NROW
      do j=1,NCOL
        if (BOTM(j,i,kbot) < this%BotMin) this%BotMin = BOTM(j,i,kbot)
      enddo
    enddo
    !
    ! If BCF is active and LAYCON=1 for layer 1, redefine TopNew.
    if (associated(LAYCON)) then
      if (LAYCON(1) == 1) then
        msg = 'The MODFLOW 6 Node Property Flow package has no equivalent' &
           // ' for LAYCON=1 of the BCF package.  As a result, calculated' &
           // ' heads may differ if heads in layer 1 are above the TOP' &
           // ' elevation.'
        call store_note(msg)
      endif
    endif
    !
    return
    end subroutine DefineTopNew

end module DisWriterModule

