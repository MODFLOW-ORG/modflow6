module FhbPackageWriterModule

  use ChdModule, only: ChdType, AddChdToList, GetChdFromList, &
                       ConstructChdType
  use ChdPackageWriterModule, only: ChdPackageWriterType
  use ConstantsModule, only: DONE, LENTIMESERIESNAME, MAXCHARLEN, &
                             DZERO
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use FileListModule, only: FileListType
  use FileTypeModule, only: FileType
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: IOUT, NCOL, NROW, NLAY, NPER, IBOUND, STRT, PERLEN
  use GlobalVariablesModule, only: masteridomain, verbose, echo
  use GWFBASMODULE, only: SGWF2BAS7PNT
  use GWFFHBMODULE, only: NBDTIM, NFLW, NHED, IFHBCB, IFLLOC, IHDLOC, &
                          BDTIM, FLWRAT, SBHED, cnstmf, cnstmh
  use GwfFhbSubs, only: GWF2FHB7AR
  use InputOutputModule, only: GetUnit
  use LineListModule, only: LineListType
  use ListModule, only: ListType
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: store_note
  use WelPackageWriterModule, only: WelPackageWriterType

  implicit none

  private
  public :: FhbPackageWriterType

  type, extends(PackageWriterType) :: FhbPackageWriterType
    ! This class is for handling FHB package stuff.
    ! It contains:
    !   WelPackageWriterType for the flow boundaries, and
    !   ChdWriterType for the head boundaries.
    !
    type (WelPackageWriterType), pointer :: WellWriter => null()
    type (ChdPackageWriterType), pointer :: ChdWriter => null()
    character(len=LENTIMESERIESNAME), allocatable, dimension(:) :: flowtsnames, &
                                                                   headtsnames
  contains
    procedure, public :: AllocatePointers
    procedure, public :: MyType
    procedure, public :: PrepareFiles
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessStressLoop
    procedure, public :: RemoveChdDuplicates
    procedure, public :: WriteFiles
    procedure, public :: WriteStressPeriodListData
  end type FhbPackageWriterType

  contains

  subroutine AllocatePointers(this)
    ! Overrides PackageWriterType%AllocatePointers
    implicit none
    ! dummy
    class(FhbPackageWriterType) :: this
    !
    allocate(this%ICbc)
    this%ICbc = 0
    allocate(this%ipr)
    allocate(this%NBndPeriod)
    allocate(this%IPB)
    allocate(this%NBNDVL)
    allocate(this%Nlaynew)
    !
    return
  end subroutine AllocatePointers

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(FhbPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: i, IFHBD3, j, k, m, nhedused, km
    character(len=MAXCHARLEN) :: fname,msg
    character(len=12) :: method
    type(WelPackageWriterType), pointer :: w => null()
    type(ChdPackageWriterType), pointer :: c => null()
    ! format
    10 format(a,'_L',i0,'_R',i0,'_C',i0,'_',i0)
    ! format
    20 format(a,' boundary removed at (',i0,',',i0,',',i0,')')
    !
    call this%AllocatePointers()
    !
    this%NAux = 0
    this%Igrid = igrid
    !
    if (this%IuOrig > 0) then
      ! Read FHB-package input file.
      call GWF2FHB7AR(this%IuOrig, igrid)
      !
      if (NFLW > 0) then
        ! Specified flows
        allocate(this%WellWriter)
        w => this%WellWriter
        w%Active = .true.
        call w%AllocatePointers()
        w%ipr = 1
        fname = trim(this%ModelBasename) // '.fhb.wel6'
        call w%InitializeFile(fname, 'WEL6', this%PackageName)
        w%PrintInput = .true.
!        w%PackageName = 'FHB'
        w%fileobj%PkgName = w%PackageName
        !
        ! Ned todo: Populate WellWriter as needed.
        w%MaxActiveBnd = NFLW
        allocate(w%NBNDVL)
        w%NBNDVL = 4
        w%NAux = 0
        allocate(w%rlist(w%NBNDVL,NFLW))
        allocate(this%flowtsnames(nflw))
        do m=1,NFLW
          k = iflloc(1,m)   ! layer
          i = iflloc(2,m)   ! row
          j = iflloc(3,m)   ! column
          if (associated(masteridomain)) then
            if (masteridomain(j,i,k)==0) then
              if (verbose) then
                msg = ''
                write(msg,20)trim(this%PkgType),k,i,j
                write(*,*)trim(msg)
              endif
              cycle
            endif
          endif
          write(this%flowtsnames(m),10)'FLOW',k,i,j,m
          w%rlist(1,m) = k
          w%rlist(2,m) = i
          w%rlist(3,m) = j
        enddo
        !
        ! Create a time-series file for flows
        w%TsFile = trim(this%ModelBasename) // '.fhb.wel6.ts'
        IFHBD3 = size(FLWRAT,1)
        if (NBDTIM == 1) then
          method = 'STEPWISE'
        else
          method = 'LINEAR'
        endif
        call write_timeseriesfile(w%TsFile, NFLW, NBDTIM, method, &
                                  this%flowtsnames, bdtim, IFHBD3, &
                                  FLWRAT, cnstmf, this%Igrid)
      endif
      !
      if (NHED > 0) then
        ! Specified heads
        c => this%ChdWriter
        c%Active = .true.
        c%source = 'FHB'
        call c%AllocatePointers()
        c%ipr = 1
        fname = trim(this%ModelBasename) // '.fhb.chd'
        call c%InitializeFile(fname, 'CHD6', this%PackageName)
        c%PrintInput = .true.
        c%PackageName = 'CHD-FHB'
        c%fileobj%PkgName = c%PackageName
        !
        ! Ned todo: populate ChdWriter as needed.
        c%MaxActiveBnd = NHED
        nhedused = NHED
        allocate(c%NBNDVL)
        c%NBNDVL = 4
        c%NAux = 0
        allocate(c%rlist(c%NBNDVL,NHED))
        allocate(c%ichdvar(NCOL,NROW,NLAY))
        c%ichdvar = 0
        allocate(this%headtsnames(NHED))
        km = 0
        do m=1,NHED
          k = IHDLOC(1,m)   ! layer
          i = IHDLOC(2,m)   ! row
          j = IHDLOC(3,m)   ! column
          c%ichdvar(j,i,k) = 1
          if (associated(masteridomain)) then
            if (masteridomain(j,i,k)==0) then
              IHDLOC(1,m) = -1
              nhedused = nhedused - 1
              if (verbose) then
                msg = ''
                write(msg,20)trim(this%PkgType),k,i,j
                write(*,*)trim(msg)
              endif
              cycle
            endif
          else
            if (associated(this%Idomain)) then
              if (this%Idomain(j,i,k)==0) then
                IHDLOC(1,m) = -1
                nhedused = nhedused - 1
                if (verbose) then
                  msg = ''
                  write(msg,20)trim(this%PkgType),k,i,j
                  write(*,*)trim(msg)
                endif
                cycle
              endif
            endif
          endif
          km = km + 1
          write(this%headtsnames(m),10)'HEAD',k,i,j,m
          c%rlist(1,m) = k
          c%rlist(2,m) = i
          c%rlist(3,m) = j
        enddo
        c%MaxActiveBnd = nhedused
        !
        ! Create a time-series file for heads
        c%TsFile = trim(this%ModelBasename) // '.fhb.chd.ts'
        IFHBD3 = size(SBHED,1)
        if (NBDTIM == 1) then
          method = 'STEPWISE'
        else
          method = 'LINEAREND'
        endif
        call write_timeseriesfile(c%TsFile, NHED, NBDTIM, method, &
                                  this%headtsnames, bdtim, IFHBD3, &
                                  SBHED, cnstmh, this%Igrid)
      endif
    endif
    !
    return
  end subroutine ProcessAllocate

  subroutine PrepareFiles(this, Mf6Files)
    implicit none
    ! dummy
    class(FhbPackageWriterType) :: this
    type(FileListType), intent(inout) :: Mf6Files
    ! local
    integer :: iu15
    character(len=MAXCHARLEN) :: fname15
    character(len=12) :: filtyp15
    !
    if (associated(this%WellWriter)) then
      if (this%WellWriter%Active) then
        fname15 = this%WellWriter%fileobj%FName
        filtyp15 = this%WellWriter%fileobj%FType
        iu15 = this%WellWriter%fileobj%IUnit
        this%WellWriter%fileobj%FCode = FCINPUT
        call Mf6Files%AddFile(fname15, filtyp15, iu15, FCINPUT, &
                                 this%WellWriter%fileobj%Pkgname)
      endif
    endif
    !
    if (associated(this%ChdWriter)) then
      if (this%ChdWriter%Active) then
        fname15 = this%ChdWriter%fileobj%FName
        filtyp15 = this%ChdWriter%fileobj%FType
        iu15 = this%ChdWriter%fileobj%IUnit
        this%ChdWriter%fileobj%FCode = FCINPUT
        call Mf6Files%AddFile(fname15, filtyp15, iu15, FCINPUT, &
                                 this%ChdWriter%fileobj%PkgName)
      endif
    endif
    !
    return
  end subroutine PrepareFiles

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(FhbPackageWriterType) :: this
    integer, intent(in) :: igrid
    !
    ! FHB does not read data by stress period, so return
    return
    !
  end subroutine ProcessStressLoop

  subroutine WriteFiles(this, igrid)
    implicit none
    ! dummy
    class(FhbPackageWriterType), intent(inout) :: this
    integer, intent(in) :: igrid
    ! formats
    20 format(2x,i0,2x,i0,2x,i0,2x,g15.8)
    60 format(a)
    70 format(/,a)
    !
    if (.not. this%Active) return
    !
    echo = .true.
    !
    if (associated(this%WellWriter)) then
      if (this%WellWriter%Active) then
        call this%WellWriter%WriteFileUsingTS(igrid, NFLW, this%flowtsnames)
      endif
    endif
    !
    if (associated(this%ChdWriter)) then
      if (this%ChdWriter%Active) then
        call this%ChdWriter%WriteFileUsingTS(igrid, NHED, this%headtsnames)
      endif
    endif
    !
    return
  end subroutine WriteFiles

  subroutine RemoveChdDuplicates(this, chdList)
    ! dummy
    class(FhbPackageWriterType) :: this
    type(ListType), pointer :: chdList
    ! local
    integer :: i1, i2, j1, j2, k, k1, k2, kc
    integer :: m, n, nchdib
    type(ListType), pointer :: templist
    type(ChdType), pointer :: chd, newchd
    character(len=100) :: msg
    ! format
    10 format(i0,' CHDs from IBOUND were removed because of duplication', &
              ' in FHB package.')
    !
    nchdib = chdList%Count()
    allocate(templist)
    k = 0
    kc = 0
    !
    ! Build new list of CHDs, keeping only non-duplicates
    outerloop: do n=1,nchdib
      chd => GetChdFromList(chdList, n)
      k2 = chd%klay
      i2 = chd%irow
      j2 = chd%jcol
      innerloop: do m=1,NHED
        k1 = IHDLOC(1,m)   ! layer
        i1 = IHDLOC(2,m)   ! row
        j1 = IHDLOC(3,m)   ! column
        if (k1 == k2 .and. i1 == i2 .and. j1 == j2) then
          k = k + 1
          cycle outerloop
        endif
      enddo innerloop
      allocate(newchd)
      call chd%CopyTo(newchd)
      kc = kc + 1
      call AddChdToList(templist, newchd)
    enddo outerloop
    !
    ! Clear and deallocate original ChdList
    call chdList%Clear(.true.)
    deallocate(chdList)
    ! Point ChdList to new list
    chdList => templist
    !
    if (k > 0) then
      write(msg,10)k
      call store_note(msg)
    endif
    !
    return
  end subroutine RemoveChdDuplicates

  subroutine WriteStressPeriodListData(this, lineList)
    ! Override if needed
    implicit none
    ! dummy
    class(FhbPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    !
    return
  end subroutine WriteStressPeriodListData

  subroutine write_timeseriesfile(tsfile, nts, ntimes, method, tsnames, &
                                  times, IFHBD3, values, cnstm, igrid)
    implicit none
    ! dummy
    character(len=*), intent(in) :: tsfile
    integer,          intent(in) :: nts, ntimes
    character(len=*), intent(in) :: method
    character(len=LENTIMESERIESNAME), dimension(nts), intent(in) :: tsnames
    double precision, dimension(ntimes),              intent(in) :: times
    integer,                                          intent(in) :: IFHBD3
    double precision, dimension(IFHBD3,nts),          intent(in) :: values
    double precision,                                 intent(in) :: cnstm
    integer, intent(in) :: igrid
    ! local
    integer :: i, iu, j
    character(len=:), allocatable :: aline
    character(len=20) :: cvalue
    double precision :: timemax
    ! formats
    1 format()
    10 format(a)
    15 format(2x,a)
    20 format(2x,a,2x,a)
    25 format(g15.8)
    30 format(2x,a,2x,g14.7)
    !
    ! Open file for writing
    iu = GetUnit()
    !
    ! Write ATTRIBUTES block
    open(iu,file=tsfile,status='REPLACE')
    write(iu,10)'BEGIN ATTRIBUTES'
    aline = 'NAME'
    do i=1,nts
      aline = aline // ' ' // trim(tsnames(i))
    enddo
    write(iu,15)aline
    write(iu,20)'METHOD', trim(method)
    if (cnstm /= DONE) then
      write(iu,30)'SFAC',cnstm
    endif
    write(iu,10)'END ATTRIBUTES'
    !
    write(iu,1)
    write(iu,10)'BEGIN TIMESERIES'
    ! Write time-series records
    do i=1,ntimes
      write(cvalue,25)times(i)
      aline = trim(cvalue)
      do j=1,nts
        write(cvalue,25)values(i,j)
        aline = aline // ' ' // trim(cvalue)
      enddo
      write(iu,10)aline
    enddo
    !
    ! Special case for a single step:
    ! Repeat line of data, but with maximum simulation time
    if (ntimes == 1 .and. method == 'STEPWISE') then
      !
      ! Determine max simulation time
      call SGWF2BAS7PNT(igrid)
      timemax = DZERO
      do i=1,NPER
        timemax = timemax + PERLEN(i)
      enddo
      !
      write(cvalue,25)timemax
      aline = trim(cvalue)
      do j=1,nts
        write(cvalue,25)values(1,j)
        aline = aline // ' ' // trim(cvalue)
      enddo
      write(iu,10)aline
    endif
    !
    write(iu,10)'END TIMESERIES'
    !
    close(iu)
    !
    return
  end subroutine write_timeseriesfile

  function MyType(this) result (ctype)
    ! dummy
    class(FhbPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'FhbPackageWriterType'
    !
    return
  end function MyType

end module FhbPackageWriterModule
