module StoWriterModule

  use ConstantsModule, only: DONE, DZERO
  use ConstantsPHMFModule, only: FCINPUT
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: ncol, nrow, nlay, nper, LAYCBD, LENUNI
  use GlobalVariablesModule, only: echo
  use SimPHMFModule, only: count_errors, store_error, store_note, ustop
  use UtilitiesModule, only: Write2Drel

  type, extends(FileWriterType) :: StoWriterType
    integer :: Istocb  ! output flows (-1, 0, 1)
    integer :: Isfac   ! indicates if ss is read as storativity
    ! Isfac = 0 by default; Isfac = 1 if STORAGECOEFFICIENT is turned on.
    integer :: Inewton ! newton-raphson flag (will not be used by converter)
    integer,               pointer :: Nlaynew => null()
    integer, dimension(:), pointer :: Layptr => null()
    integer, pointer, dimension(:)  :: Iconvert  ! confined (0) or convertible (1), upstream (4)?
    double precision, pointer, dimension(:,:,:) :: Ss => null() ! specific storage (or storage coefficient)
    double precision, pointer, dimension(:,:,:) :: Sy => null() ! specific yield
    logical, pointer, dimension(:)  :: Transient => null()
    logical :: Newton
  contains
    procedure :: AllocateArrays
    procedure :: InitializeFile
    procedure :: WriteFile
    procedure :: WriteOptions
    procedure :: WriteStoData
  end type StoWriterType

contains

  subroutine AllocateArrays(this)
    implicit none
    class(StoWriterType) :: this
    !
    allocate(this%Iconvert(nlay))
    allocate(this%Ss(ncol,nrow,nlay))
    allocate(this%Sy(ncol,nrow,nlay))
    allocate(this%Transient(nper))
    this%Ss = 0.0
    this%Sy = 0.0
    this%Transient = .false.
    this%Newton = .false.
    !
    return
  end subroutine AllocateArrays

  subroutine InitializeFile(this, fname, ftype, pkgname)
    implicit none
    ! dummy
    class(StoWriterType), intent(inout) :: this
    character(len=*), intent(in)    :: fname
    character(len=*), intent(in)    :: ftype
    character(len=*), intent(in), optional :: pkgname
    ! local
    !
    ! Invoke superclass initializer
    if (present(pkgname)) then
      call this%FileWriterType%InitializeFile(fname, ftype, pkgname)
    else
      call this%FileWriterType%InitializeFile(fname, ftype)
    endif
    !
    ! Assign file code
    this%fileobj%FCode = FCINPUT
    !
    ! Assign defaults
    this%Isfac = 0 ! Ss is read (Isfac=1 means S is read)
    !
    return
  end subroutine InitializeFile

  subroutine WriteFile(this)
    implicit none
    ! dummy
    class(StoWriterType) :: this
    ! local
    logical :: iconvdone, trans
    integer :: iprni, iprnr, iu, n
    ! formats
    5  format()
    54 format(4x,a)
    70 format('BEGIN Period',1x,i0)
    72 format('END Period')
    !
    if (.not. this%Active) return
    if (echo) then
      iprni = 2
      iprnr = 12
    else
      iprni = -2
      iprnr = -12
    endif
    !
    call this%WriteOptions()
    call this%WriteStoData(iprnr)
    iu = this%fileobj%IUnit
    iconvdone = .false.
    !
    ! Assume transient, as STO package does.
    trans = .true.
    !
    do n=1,nper
      if (trans .neqv. this%Transient(n)) then
        ! Switch from steady to transient or vice versa this period
        trans = this%Transient(n)
        ! Write BEGIN Period
        write(iu,5)
        write(iu,70)n
        ! Write TRANSIENT or STEADY-STATE
        if (this%Transient(n)) then
          write(iu,54)'TRANSIENT'
        else
          write(iu,54)'STEADY-STATE'
        endif
        ! Write END Period
        write(iu,72)
      endif
    enddo
    !
    return
  end subroutine WriteFile

  subroutine WriteOptions(this)
    implicit none
    ! dummy
    class(StoWriterType) :: this
    ! local
    integer :: iu
    ! formats
    5  format()
    10 format('BEGIN Options')
    50 format(2x,a)
    100 format('END Options')
    !
    iu = this%fileobj%IUnit
    !
    ! Write BEGIN Options
    write(iu,5)
    write(iu,10)
    !
    if (this%Istocb<0) write(iu,50)'SAVE_FLOWS'
!    if (this%Istocb>0) write(iu,50)'PRINT_FLOWS' ! not yet supported??
    if (this%Isfac/=0) write(iu,50)'STORAGECOEFFICIENT'
!    if (.not. this%Newton) write(iu,50)'NO_NEWTON'
    !
    ! Write END Options
    write(iu,100)
    !
    return
  end subroutine WriteOptions

  subroutine WriteStoData(this, iprnr)
    implicit none
    ! dummy
    class(StoWriterType) :: this
    integer, intent(in) :: iprnr
    ! local
    integer :: i, iu, j, k, knew
    logical :: constant
    character(len=500) :: msg
    double precision :: rconst, ss, sy
    ! formats
    5  format()
    10 format(a)
    30 format(4x,a,2x,i0)
    40 format(4x,a,2x,g15.8)
    43 format(4x,a,2x,g10.3)
    50 format(2x,a)
    100 format(a,1x,i0,1x,a,1x,i0,1x,a,1x,g10.3)
    !
    ! Write BEGIN StoData
    iu = this%fileobj%IUnit
    write(iu,5)
    write(iu,10)'BEGIN GRIDDATA'
    !
    ! Iconvert
    write(iu,50)'ICONVERT LAYERED'
    do k=1,nlay
      write(iu,30)'CONSTANT', this%Iconvert(k)
      if (LAYCBD(k) /= 0) then
        ! insert an entry for a quasi-3d unit -- assume nonconvertible
        write(iu,30)'CONSTANT', 0
      endif
    enddo
    !
    ! Ss (or S)
    write(iu,50)'SS LAYERED'
    knew = 0
    do k=1,nlay
      knew = knew + 1
      constant = .true.
      rconst = this%Ss(1,1,k)
      checkss: do i=1,nrow
        do j=1,ncol
          if (this%Ss(j,i,k)/=rconst) then
            constant = .false.
            exit checkss
          endif
        enddo
      enddo checkss
      if (constant) then
        write(iu,40)'CONSTANT', rconst
      else
        call Write2Drel(iu, nrow, ncol, this%Ss(:,:,k), &
                        constant, rconst, 'SS', .false., iprnr)
      endif
      if (LAYCBD(k) /= 0) then
        ! insert an entry for a quasi-3d unit.
        ss = DZERO
!              ss = get_q3d_ss()
        write(iu,43)'CONSTANT', ss
        knew = knew + 1
        write(msg,100)'Original quasi-3D unit below layer',k, &
               'is now layer',knew,'and has been assigned SS value:',ss
       call store_note(msg)
      endif
    enddo
    !
    ! Sy
    write(iu,50)'SY LAYERED'
    knew = 0
    do k=1,nlay
      knew = knew + 1
      constant = .true.
      rconst = this%Sy(1,1,k)
      checksy: do i=1,nrow
        do j=1,ncol
          if (this%Sy(j,i,k)/=rconst) then
            constant = .false.
            exit checksy
          endif
        enddo
      enddo checksy
      if (constant) then
        write(iu,40)'CONSTANT', rconst
      else
        call Write2Drel(iu, nrow, ncol, this%Sy(:,:,k), &
                        constant, rconst, 'SY', .false., iprnr)
      endif
      if (LAYCBD(k) /= 0) then
        ! insert an entry for a quasi-3d unit.
!              sy = get_q3d_sy()
        sy = DZERO
        write(iu,43)'CONSTANT', sy
        knew = knew + 1
        write(msg,100)'Original quasi-3D unit below layer',k, &
               'is now layer',knew,'and has been assigned SY value:',sy
       call store_note(msg)
      endif
    enddo
    !
    ! Write END StoData
    write(iu,10)'END GRIDDATA'
    !
    return
  end subroutine WriteStoData

  double precision function get_q3d_ss() result(ss)
    ! Based on LENUNI, return a reasonable, small value for
    ! specific storage, to be used in active cells that will
    ! replace a quasi-3D confining unit.
    implicit none
    ! local
    ! -- ssms is a typical, small specific storage in 1/meter.
    double precision :: converter, ssms
    character(len=300) :: msg
    character(len=10)  :: strng
    ! format
    10 format(g10.3)
    !
    ! -- Start with a small Ss value (1/meter)
    ssms = 1.0e-9
    !
    ! -- Convert to model units specified in BAS input
    if (LENUNI < 1 .or. LENUNI > 3) then
      msg = 'Error: LENUNI is undefined in Discretization Package input. ' // &
          'A valid value needs to be assigned so that SS can be ' // &
          'assigned to cells of former quasi-3D unit.'
      call store_error(msg)
    endif
    if (count_errors() > 0) call ustop()
    !
    !
    ! -- Define length conversion
    select case (LENUNI)
    case (1)
      ! 1/foot
      converter = 0.3048d0
    case (2)
      ! 1/meter -- no length conversion needed
      converter = DONE
    case (3)
      ! 1/centimeter
      converter = 1.0d-2
    end select
    !
    ! -- Perform conversion
    ss = ssms * converter
    ! -- Round the value to about 3 significant figures
    write(strng,10)ss
    read(strng,*)ss
    !
    return
  end function get_q3d_ss

  double precision function get_q3d_sy() result(sy)
    ! Return a reasonable, small value for specific yield,
    ! to be used in active cells that will replace a
    ! quasi-3D confining unit.
    implicit none
    !
    sy = 0.0001d0   ! dimensionless
    !
    return
  end function get_q3d_sy

end module StoWriterModule
