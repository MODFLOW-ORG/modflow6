module TdisWriterModule

  use ConstantsModule, only: DZERO
  use FileWriterModule, only: FileWriterType
  use SimPHMFModule, only: store_error, ustop
  use StressPeriodModule, only: StressPeriodType

  implicit none

  private
  public :: TdisWriterType

  type, extends(FileWriterType) :: TdisWriterType
    character(len=10), private :: TimeUnit = 'UNDEFINED'
    integer, private :: Nper = -1
    type(StressPeriodType), dimension(:), pointer :: StressPeriods => null()
  contains
    procedure, public :: InitializeFile => initialize_tdis
    procedure, public :: Alloc
    procedure, public :: WriteFile
    procedure, public :: AssignTimeUnit
    procedure, public :: GetStressPeriodTimes
  end type TdisWriterType

contains

  subroutine initialize_tdis(this, fname, ftype, pkgname)
    implicit none
    ! dummy
    class(TdisWriterType), intent(inout) :: this
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
    ! IC is always active
    this%Active = .true.
    !
    return
  end subroutine initialize_tdis

  subroutine Alloc(this, nper)
    implicit none
    ! dummy
    class(TdisWriterType), intent(inout) :: this
    integer, intent(in) :: nper
    ! local
    character(len=100) :: msg
    !
    if (nper > 0) then
      this%Nper = nper
      if (.not. associated(this%StressPeriods)) then
        allocate(this%StressPeriods(nper))
      endif
    else
      msg = 'Error in TdisWriterType%Alloc: NPER <= 0'
      call store_error(msg)
      call ustop()
    endif
    return
  end subroutine Alloc

  subroutine WriteFile(this)
    implicit none
    ! dummy
    class(TdisWriterType), intent(inout) :: this
    ! local
    integer :: i, iu
    ! formats
    5  format()
    10 format('BEGIN ',a)
    20 format(2x,'TIME_UNITS ',a)
    30 format(2x,'NPER ',i0)
    40 format(2x,g15.8,2x,i10,2x,g15.8,2x,'Items: PERLEN NSTP TSMULT')
    50 format('END ',a)
    !
    iu = this%fileobj%IUnit
    ! Write Options block
    write(iu,5)
    write(iu,10)'Options'
    write(iu,20)trim(this%TimeUnit)
    write(iu,50)'Options'
    write(iu,5)
    !
    ! Write Dimensions block
    write(iu,10)'Dimensions'
    write(iu,30)this%nper
    write(iu,50)'Dimensions'
    write(iu,5)
    !
    ! Write Stress_Periods block
    write(iu,10) 'PERIODDATA'
    do i=1,this%Nper
      write(iu,40)this%StressPeriods(i)%perlen, this%StressPeriods(i)%nstp, &
                  this%StressPeriods(i)%tsmult
    enddo
    write(iu,50) 'PERIODDATA'
    !
    return
  end subroutine WriteFile

  subroutine AssignTimeUnit(this, TimeUnit)
    implicit none
    class(TdisWriterType), intent(inout) :: this
    character(len=*), intent(in) :: TimeUnit
    !
    this%TimeUnit = TimeUnit
    return
  end subroutine AssignTimeUnit

  subroutine GetStressPeriodTimes(this, kper, starttime, endtime)
    ! dummy
    class(TdisWriterType) :: this
    integer, intent(in) :: kper
    double precision, intent(out) :: starttime, endtime
    ! local
    integer :: iper
    !
    starttime = DZERO
    do iper = 1,kper-1
      starttime = starttime + this%StressPeriods(iper)%perlen
    enddo
    endtime = starttime + this%StressPeriods(kper)%perlen
    !
    return
  end subroutine GetStressPeriodTimes

end module TdisWriterModule
