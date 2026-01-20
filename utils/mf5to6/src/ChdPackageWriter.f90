module ChdPackageWriterModule

  use ChdModule, only: ChdType, CastAsChdType, ConstructChdType, &
                       AddChdToList, GetChdFromList
  use ChdObsWriterModule, only: ChdObsWriterType, createChdObsWriter
  use ConstantsModule, only: MAXCHARLEN, LINELENGTH, &
                             DZERO, LENTIMESERIESNAME, LENHUGELINE, DONE, &
                             DEM1, DEM6, DNODATA
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use FileTypeModule, only: FileType
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: IOUT, NCOL, NROW, NLAY, NPER, IBOUND, STRT, TSMULT, NSTP, &
                    PERLEN
  use GlobalVariablesModule, only: masteridomain, verbose, echo
  use GWFCHDMODULE, only: NCHDS, MXCHD, NCHDVL, IPRCHD, NPCHD, &
                          ICHDPB, NNPCHD, CHDAUX, CHDS, GWF2CHD7DA
  use GwfChdSubs, only: GWF2CHD7AR, GWF2CHD7RP
  use GWFFHBMODULE, only: IHDLOCfhb => IHDLOC, NHEDhfb => NHED
  use InputOutputModule, only: GetUnit, openfile, URWORD
  use LineListModule, only: LineListType
  use ListModule, only: ListType
  use PackageWriterModule, only: PackageWriterType
  use SimPHMFModule, only: store_note, store_error, ustop
  use TdisVariablesModule, only: GlobalTdisWriter
  use TimeSeriesRecordModule, only: TimeSeriesRecordType, &
                                    ConstructTimeSeriesRecord
  use TimeSeriesModule, only: TimeSeriesFileType, SameTimeSeries
  use UtilitiesModule, only: close_file

  implicit none

  private
  public :: ChdPackageWriterType

  type, extends(PackageWriterType) :: ChdPackageWriterType
    ! This class is for handling constant heads defined by
    ! negative values in the IBOUND array. It also
    ! reads and processes CHDs defined by stress period
    ! in a CHD-package input file.
    ! ichdvar <0: constant by stress period; >0: varying, need time series
    integer :: kpercurrent
    integer :: numtimeseries = 0
    integer :: maxchds = 0
    integer,          dimension(:,:,:), pointer :: ibnd => null()
    integer,          dimension(:,:,:), pointer :: ichdvar => null()
    double precision, dimension(:,:,:), pointer :: chdshead => null()
    double precision, dimension(:,:,:), pointer :: chdehead => null()
    type(ListType),            pointer :: IbChdList => null()                    ! CHDs from IBOUND
    type(ListType),            pointer :: TvChdList => null()                    ! Time-variable CHDs from CHD package
    type(ListType),            pointer :: CnstChdList => null()                  ! Constant CHDs from CHD package
    type(TimeSeriesFileType), pointer :: ChdTsGroup => null()
    !type(ChdObsWriterType),    pointer :: ChdObsWriter => null()
  contains
    ! Public procedures
    procedure, public :: AllocatePointers
    procedure, public :: InitializeFile
    procedure, public :: MyType
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessIbound
    procedure, public :: ProcessStressLoop
    procedure, public :: WriteDimensions
    procedure, public :: WriteFile
    procedure, public :: WriteFileUsingTS
    procedure, public :: WriteStressPeriodListData
    ! Private procedures
    procedure, private :: ExtractTimeSeriesData
    procedure, private :: GetMaxbound
    procedure, private :: InflateTimeSeries
    procedure, private :: PopulateIchdvar
    procedure, private :: ReduceTimeSeries
    procedure, private :: RewriteOutputFile
    procedure, private :: SubstituteTsName
    procedure, private :: WriteConstantListData
    procedure, private :: WriteIboundChdListData
    procedure, private :: WriteTimeSeriesFile
    procedure, private :: WriteTimeVaryingListData
  end type ChdPackageWriterType

contains

  ! Type-bound public procedures

  subroutine AllocatePointers(this)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    ! local
    type(ChdObsWriterType), pointer :: newChdObsWriter => null()
    !
    ! Call parent class method
    call this%PackageWriterType%AllocatePointers()
    !
    ! Allocate lists of CHD cells and ObsWriter member
    allocate(this%IbChdList)
    allocate(this%TvChdList)
    allocate(this%CnstChdList)
    call createChdObsWriter(newChdObsWriter, this%ModelBasename, this%IuOrig)
    if (associated(this%PkgObsWriter)) then
      deallocate(this%PkgObsWriter)
    endif
    newChdObsWriter%PkgIbChdList => this%IbChdList
    newChdObsWriter%PkgTvChdList => this%TvChdList
    newChdObsWriter%PkgCnstChdList => this%CnstChdList
    this%PkgObsWriter => newChdObsWriter
    !
    return
  end subroutine AllocatePointers

  subroutine InitializeFile(this, fname, ftype, pkgname)
    implicit none
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    character(len=*), intent(in) :: fname
    character(len=*), intent(in) :: ftype
    character(len=*), intent(in), optional :: pkgname
    !
    ! Invoke superclass initializer
    if (present(pkgname)) then
      call this%FileWriterType%InitializeFile(fname, ftype, pkgname)
    else
      call this%FileWriterType%InitializeFile(fname, ftype, this%PackageName)
    endif
    !
    ! Assign pointer
    this%ibnd => IBOUND
    !
    ! Initialize options
    this%PrintInput = .true.
    this%PrintFlows = .true.
    this%SaveFlows = .false.
    this%PackageName = ''
    this%PkgType = 'CHD'
    this%DefaultBudgetText = 'CONSTANT HEAD'
    !
    return
  end subroutine InitializeFile

  function MyType(this) result (ctype)
    ! dummy
    class(ChdPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'ChdPackageWriterType'
    !
    return
  end function MyType

  subroutine ProcessAllocate(this, igrid)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: fname, msg
    !
    !call this%AllocatePointers()
    allocate(this%ChdTsGroup)
    this%ipr = 1

    this%Active = .true.
    this%DefaultBudgetText = 'CONSTANT HEAD'
    ! Ned todo: not sure that this will work for LGR--need to look at LGR source
    this%fileobj%FCode = FCINPUT
    !
    if (this%source == 'CHD') then
      this%fileobj%FType = 'CHD6'
      this%PkgType = 'CHD'
      fname = trim(this%ModelBasename) // '.chd'
      call this%FileWriterType%InitializeFile(fname, this%fileobj%FType, &
                                              this%PackageName)
    endif

    this%NStressDim = 2   ! Specific to CHD package; for Shead, Ehead.
    this%NAux = 0
    !
    if (this%IuOrig > 0) then
      ! Allocate ichdvar, populate ichdvar by reading entire CHD input file,
      ! then rewind CHD file.
      allocate(this%ichdvar(NCOL,NROW,NLAY))
      call this%PopulateIchdvar(igrid)
      ! Dimension chdshead and chdehead according to layering original in
      ! original model.
      allocate(this%chdshead(NCOL,NROW,NLAY))
      allocate(this%chdehead(NCOL,NROW,NLAY))
      this%chdshead = DNODATA
      this%chdehead = DNODATA
      !
      ! Read CHD-package input file (but parameter data have already been
      ! stored, so don't store again).
      call GWF2CHD7AR(this%IuOrig, igrid, .false.)
      !
      ! NStressDim is number of values after k,i,j required to
      ! specify boundary stress.
      !
      ! NVALP1 is rlist index for element following minimum required data.
      this%nvalp1 = this%NStressDim + 4
      !
      ! Assign package-specific pointers
      this%rlist => CHDS
      this%ipr => IPRCHD
      allocate(this%ICbc)
      this%ICbc = 0
      this%NBndPeriod => NCHDS
      this%Aux => CHDAUX
      this%IPB => ICHDPB
      this%MaxActiveBnd = this%IPB - 1
      this%NBNDVL => NCHDVL
      this%NAux = this%NBNDVL - this%nvalp1 + 1
      if (this%NAux > 0) then
        msg = 'Auxiliary variables are not supported in conversion of CHD package.'
        call store_note(msg)
        this%NAux = 0
      endif
      this%nstop = this%nvalp1 - 1 + this%NAux
    endif
    !
    ! define format for printing boundary data
!    write(ctemp,'(i0)')this%NAux + this%NStressDim
!    this%fmat = '(3(2x,i0),' // trim(ctemp) // '(2x,g15.8))'
     this%fmat = '(3(2x,i0),2x,a)'
    !
    if (this%source == 'CHD') then
      call this%WriteOptions()
      call this%WriteDimensions()
    endif
    !
    return
  end subroutine ProcessAllocate

  subroutine ProcessIbound(this)
    implicit none
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    ! local
    integer :: i, j, k, kk, kold, n
    type(ChdType), pointer :: chd => null()
    class(*),      pointer :: obj => null()
    character(len=MAXCHARLEN) :: fname
    ! formats
    10 format(/,i0, &
        ' CHD nodes have been defined from IBOUND and starting heads')
    !
    kk = 0
    fname = trim(this%ModelBasename) // '.ibound.chd'
    ! Iterate through IBOUND to find CH cells (IBOUND < 0)
    do kold=1,NLAY
      k = this%Layptr(kold)
      do i=1,NROW
        do j=1,NCOL
          if (this%ibnd(j,i,kold)<0) then
            if (.not. this%Active) then
              this%Active = .true.
              call this%InitializeFile(fname, 'CHD6', this%PackageName)
              this%fileobj%FCode = FCINPUT
              this%PackageName = 'CHD-FROM-IBOUND'
              this%fileobj%PkgName = this%PackageName
            endif
            ! Create a CH node
            call ConstructChdType(chd)
            chd%irow = i
            chd%jcol = j
            chd%klay = k
            ! Define head STRT
            chd%head = STRT(j,i,kold)
            ! Add this chd to IbChdList
            obj => chd
            kk = kk + 1
            call this%IbChdList%Add(obj)
          endif
        enddo
      enddo
    enddo
    !
    n = this%IbChdList%Count()
    write(IOUT,10)n
    !
    return
  end subroutine ProcessIbound

  subroutine ProcessStressLoop(this, igrid)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: k, kper
    logical :: currentA
    character(len=50) :: procmsg
    ! formats
    10 format(/,'Processing CHD package for stress period ',i0)
    !
    ! Initially, current block is BlockA; alternate each stress period
    this%CurrentBlock => this%BlockA
    this%PreviousBlock => this%BlockB
    currentA = .true.
    procmsg = 'Processing CHD package data from ' // trim(this%source) // ' input...'
    !
    do kper=1,nper
      if (kper==1) write(*,*)trim(procmsg)
      write(iout,10)kper
      ! Read MF2005 input for current stress period
      if (this%IuOrig > 0) then
        call GWF2CHD7RP(this%IuOrig, igrid)
        call this%ExtractTimeSeriesData(kper)
      endif
      ! Write MF6 input for current stress period to a LineList
      call this%CurrentBlock%Clear(.true.)
      this%kpercurrent = kper
      call this%WriteStressPeriodListData(this%CurrentBlock)
      ! Write block to MF6 input file if needed
      call this%WriteBlockIfNeeded(kper)
      !
      ! switcheroo
      if (currentA) then
        this%CurrentBlock => this%BlockB
        this%PreviousBlock => this%BlockA
        currentA = .false.
      else
        this%CurrentBlock => this%BlockA
        this%PreviousBlock => this%BlockB
        currentA = .true.
      endif
    enddo
    !
    k = this%TvChdList%Count()
    if (k > 0) then
      ! Write all time-series files required for this ChdPackage.
      call this%WriteTimeSeriesFile()
      ! Rewrite output file to insert time-series file names into Options block.
      call this%RewriteOutputFile()
    endif
    !
    call this%BlockA%Clear(.true.)
    call this%BlockB%Clear(.true.)
    !
    return
  end subroutine ProcessStressLoop

  subroutine WriteDimensions(this)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    ! local
    integer :: i, iu, j, k, n, maxb
    ! formats
    30 format('')
    50 format(2x,a,2x,i0)
    60 format(a)
    !
    iu = this%fileobj%IUnit
    !
    maxb = this%GetMaxbound()
    !
    write(iu,30)
    write(iu,60)'BEGIN Dimensions'
    write(iu,50)'MAXBOUND', maxb
    write(iu,60)'END Dimensions'
    !
    return
  end subroutine WriteDimensions

  subroutine WriteFile(this, igrid)
    implicit none
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    integer, intent(in) :: igrid
    ! local
    ! formats
    20 format(2x,i0,2x,i0,2x,i0,2x,g15.8)
    60 format(a)
    70 format(/,a)
    !
    if (.not. this%Active) return
    !
    call this%ProcessAllocate(igrid)
    call this%WriteOptions()
    call this%WriteDimensions()
    call this%ProcessStressLoop(igrid)
    !
    call this%CloseFile()
    !
    return
  end subroutine WriteFile

  subroutine WriteFileUsingTS(this, igrid, ndim, tsnames)
    implicit none
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    integer, intent(in) :: igrid
    integer, intent(in) :: ndim
    character(len=*), dimension(ndim) :: tsnames
    ! local
    integer :: i, iu, j, k, m
    ! formats
    1 format()
    5 format(a)
    10 format(2x,3(i0,2x),a)
    !
    if (.not. this%Active) return
    !
    call this%WriteOptions()
    call this%WriteDimensions()
    !
    ! Ned todo: write CHD info here (needed for FHB-based constant heads)
    if (ndim > 0) then
      iu = this%fileobj%IUnit
      write(iu,1)
      write(iu,5)'BEGIN PERIOD 1'
      do m=1, ndim
        k = IHDLOCfhb(1,m)   ! layer
        if (k > 0) then
          i = IHDLOCfhb(2,m)   ! row
          j = IHDLOCfhb(3,m)   ! column
          write(iu,10)k, i, j, trim(tsnames(m))
        endif
      enddo
      write(iu,5)'END PERIOD'
    endif
    !
    return
  end subroutine WriteFileUsingTS

  subroutine WriteStressPeriodListData(this, lineList)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    !
    call this%WriteIboundChdListData(lineList)
    call this%WriteTimeVaryingListData(lineList)
    call this%WriteConstantListData(lineList)
    !
    return
  end subroutine WriteStressPeriodListData

  ! Type-bound private procedures

  subroutine ExtractTimeSeriesData(this, kper)
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    integer, intent(in) :: kper
    ! local
    integer, save :: kperold = 0
    integer :: i, ichd, ii, j, k, kold, nchd
    double precision :: ehead, shead, spbegintime, spendtime, &
                        spbeginplus, sheadplus, epsi, plen, tplus
    double precision, save :: totim
    double precision, parameter :: fract = 1.0d-5
    logical :: found, timevarying
    character(len=LENTIMESERIESNAME) :: tsname
    type(ChdType), pointer :: chd => null()
    type(TimeSeriesRecordType), pointer :: tsr => null()
    ! format
    10 format('CHTS_L',i0,'_R',i0,'_C',i0)
    !
    ! Find simulation times of start and end of current stress period
    call GlobalTdisWriter%GetStressPeriodTimes(kper, spbegintime, spendtime)
    plen = spendtime-spbegintime
    tplus = fract*(plen)
    spbeginplus = spbegintime + tplus
    !
    ! Calculate epsilon appropriate for differentiating
    ! between spbegintime and spbeginplus, and update totim.
    epsi = DEM6
    if (kper == 1) then
      totim = DZERO
    else
      if (totim > DZERO) then
        epsi = DEM1*(tplus/totim)
      endif
      if (epsi > DEM6) epsi = DEM6
    endif
    ! After assigning epsi, advance totim to be
    ! time at end of current stress period.
    totim = totim + plen
    !
    ! Save ehead from previous stress period as shead for current stress period
    if (kper > 1) then
      do k=1,NLAY
        do i=1,NROW
          do j=1,NCOL
            this%chdshead(j,i,k) = this%chdehead(j,i,k)
          enddo
        enddo
      enddo
    endif
    !
    ! Iterate through data read from CHD package this stress period and zero
    ! shead and ehead to support additive heads.
    do ii=1,this%NBndPeriod
      kold = nint(this%rlist(1,ii))
      i = nint(this%rlist(2,ii))  ! Row
      j = nint(this%rlist(3,ii))  ! Column
      this%chdshead(j,i,kold) = DZERO
      this%chdehead(j,i,kold) = DZERO
    enddo
    !
    ! Iterate through data read from CHD package input for current stress period.
    ! Use current data to update chdshead and chdehead for listed entries.
    !
    do ii=1,this%NBndPeriod
      kold = nint(this%rlist(1,ii))
      i = nint(this%rlist(2,ii))  ! Row
      j = nint(this%rlist(3,ii))  ! Column
      shead = this%rlist(4,ii)
      ehead = this%rlist(5,ii)
      ! Specified heads are additive
      this%chdshead(j,i,kold) = this%chdshead(j,i,kold) + shead
      this%chdehead(j,i,kold) = this%chdehead(j,i,kold) + ehead
    enddo
    !
    layers: do kold=1,NLAY
      k = this%Layptr(kold)
      rows: do i=1,NROW
        columns: do j=1,NCOL
          if (this%ichdvar(j,i,kold) /= 0) then
            ! Cell is constant head
            timevarying = (this%ichdvar(j,i,kold) > 0)
            !
            if (timevarying) then
              ! Construct/check for time series for CHD boundaries that
              ! are time-varying within at least one stress period.
              if (kper >= this%ichdvar(j,i,kold)) then
                ! CHD for this cell is active this stress period
                shead = this%chdshead(j,i,kold)
                ehead = this%chdehead(j,i,kold)
                sheadplus = shead + fract*(ehead-shead)
                !
                ! Look for existing CHD item in TvChdList for this cell
                found = .false.
                chd => null()
                nchd = this%TvChdList%Count()
                do ichd=1,nchd
                  chd => GetChdFromList(this%TvChdList, ichd)
                  if (chd%irow==i .and. chd%jcol==j .and. chd%klay==k) then
                    found = .true.
                    chd%ActiveByStressPeriod(kper:nper) = .true.
                    chd%ActiveInMf2005(kper) = .true.
                    exit
                  endif
                enddo
                !
                if (.not. found) then
                  ! Create a new time-varying CHD for this cell and add it to TvChdList
                  call ConstructChdType(chd)
                  allocate(chd%ActiveByStressPeriod(NPER))
                  chd%ActiveByStressPeriod = .false.
                  ! Once a cell is made constant-head, it stays constant-head for
                  ! rest of simulation.
                  chd%ActiveByStressPeriod(kper:nper) = .true.
                  ! Store activeness as coded in MF2005 input.
                  allocate(chd%ActiveInMf2005(NPER))
                  chd%ActiveInMf2005 = .false.
                  chd%ActiveInMf2005(kper) = .true.
                  chd%irow = i
                  chd%jcol = j
                  chd%klay = k
                  write(tsname,10)k,i,j
                  call chd%timeSeries%InitializeTimeSeries(this%ChdTsGroup,tsname)
                  chd%OriginalTsName = tsname
                  call AddChdToList(this%TvChdList, chd)
                  ! Add a time-series record for beginning of
                  ! stress period kper to the chd item.
                  call ConstructTimeSeriesRecord(tsr, spbegintime, shead)
                  call chd%timeSeries%AddTimeSeriesRecord(tsr)
                  ! Add a time-series record for a time slightly later than the
                  ! beginning of stress period kper to the chd item.
                  call ConstructTimeSeriesRecord(tsr, spbeginplus, sheadplus)
                  call chd%timeSeries%AddTimeSeriesRecord(tsr)
                  ! Add a time-series record for end of
                  ! stress period kper to the chd item.
                  call ConstructTimeSeriesRecord(tsr, spendtime, ehead)
                  call chd%timeSeries%AddTimeSeriesRecord(tsr)
                else
                  ! Add time-series records, times, and/or head values for
                  ! beginning and end of stress period to existing CHD item as needed.
                  !
                  tsr => chd%timeSeries%GetTimeSeriesRecord(spbegintime, epsi)
                  if (associated(tsr)) then
                    ! Chd item contains a time-series record for the beginning of the
                    ! stress period. Add the shead value to the existing value, but
                    ! only if it's stress period 1 or if this chd was not active
                    ! in previous stress period.
                    if (kper == 1) then
                      tsr%tsrValue = tsr%tsrValue + shead
                    elseif (.not. chd%ActiveInMf2005(kper-1)) then
                      tsr%tsrValue = tsr%tsrValue + shead
                    endif
                  else
                    ! Chd item does not contain a time-series record for the beginning
                    ! of the stress period. Make a new time-series record for the
                    ! beginning of the stress period and add it to the time series.
                    ! These instructions probably should never be executed.
                    call ConstructTimeSeriesRecord(tsr, spbegintime, shead)
                    call chd%timeSeries%AddTimeSeriesRecord(tsr)
                  endif
                  !
                  tsr => chd%timeSeries%GetTimeSeriesRecord(spbeginplus, epsi)
                  if (associated(tsr)) then
                    ! Chd item contains a time-series record for a time slightly later
                    ! than the beginning of the stress period.
                    ! Add the sheadplus value to the existing value.
                    tsr%tsrValue =  tsr%tsrValue + sheadplus
                  else
                    ! Chd item does not contain a time-series record for the time
                    ! slightly later than the beginning of the stress period.
                    ! Make a new time-series record for this time and add it
                    ! to the time series.
                    call ConstructTimeSeriesRecord(tsr, spbeginplus, sheadplus)
                    call chd%timeSeries%AddTimeSeriesRecord(tsr)
                  endif
                  !
                  tsr => chd%timeSeries%GetTimeSeriesRecord(spendtime, epsi)
                  if (associated(tsr)) then
                    ! Chd item contains a time-series record for the end
                    ! of the stress period. Add the ehead value to the existing value.
                    tsr%tsrValue =  tsr%tsrValue + ehead
                  else
                    ! Chd item does not contain a time-series record for the end
                    ! of the stress period. Make a new time-series record for the
                    ! end of the stress period and add it to the time series.
                    call ConstructTimeSeriesRecord(tsr, spendtime, ehead)
                    call chd%timeSeries%AddTimeSeriesRecord(tsr)
                  endif
                  !
                endif
              endif
            else
              ! CHD is not time-varying; add it to list of constant
              ! CHDs from CHD package.
              if (this%ichdvar(j,i,kold) < 0) then
                call ConstructChdType(chd)
                chd%klay = k
                chd%irow = i
                chd%jcol = j
                call AddChdToList(this%CnstChdList, chd)
              endif
            endif
          endif
        enddo columns
      enddo rows
    enddo layers
    !
    kperold = kper
    !
    return
  end subroutine ExtractTimeSeriesData

  function GetMaxbound(this) result (maxb)
    ! Get MAXBOUND as number of CHDs from IBOUND
    ! plus CHDs from CHD package.
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    ! local
    integer :: i, j, k, maxb
    !
    !
    ! IBOUND-based CHDs are in this%IbChdList.
    maxb = this%IbChdList%Count()
    !
    ! CHD-package based CHDs are coded in this%ichdvar
    if (associated(this%ichdvar)) then
      do k=1,NLAY
        do i=1,NROW
          do j=1,NCOL
            if (this%ichdvar(j,i,k) /= 0) then
              maxb = maxb + 1
            endif
          enddo
        enddo
      enddo
    endif
    !
    return
  end function GetMaxbound

  subroutine InflateTimeSeries(this)
    ! Ensure that time series for all time-varying CHDs include TS records
    ! for the same set of simulation times. If any records are
    ! missing, insert records with NODATA as tsrValue.
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    ! local
    integer :: i, ismall, j, nts
    integer :: kdone, lastCount
    double precision :: badtime, endtime, smallest, time
    double precision, allocatable, dimension(:) :: times, values
    logical :: allequal, lasttime
    !logical, allocatable, dimension(:) :: complete
    type(TimeSeriesRecordType), pointer :: tsr, tsRecEarlier, tsRecLater, tsrNew
    type(ChdType), pointer :: chd
    class(*), pointer :: obj => null()
    !
    lasttime = .false.
    badtime = huge(badtime)
    nts = this%TvChdList%Count()
    if (nts < 2) return
    !
    allocate(times(nts))
    allocate(values(nts))
    times = badtime
    values = DNODATA
    time = DZERO
    endtime = -1.0d0
    !
    ! Find latest time stored in all time series
    do j=1,nts
      chd => GetChdFromList(this%TvChdList, j)
      time = chd%timeSeries%FindLatestTime()
      if (time > endtime) endtime = time
    enddo
    !
    ! Reset all time series and store time and value for first node of each time series
    lastCount = 0
    do j=1,nts
      chd => GetChdFromList(this%TvChdList, j)
      call chd%timeSeries%Reset()
      tsr => chd%timeSeries%GetNextTimeSeriesRecord()
      if (associated(tsr)) then
        times(j) = tsr%tsrTime
        values(j) = tsr%tsrValue
      else
        times(j) = badtime
        values(j) = DNODATA
        lastCount = lastCount + 1
      endif
    enddo
    call get_smallest(nts, times, ismall, smallest, allequal)
    !
    ! Ensure all time series have records for current smallest time
    if (lastCount /= nts) then
      loop1: do
        if (.not. allequal) then
          do j=1,nts
            if (times(j) > smallest) then
              chd => GetChdFromList(this%TvChdList, j)
              ! Insert a NODATA TS record
              call ConstructTimeSeriesRecord(tsrNew, smallest, DNODATA)
              call chd%timeSeries%InsertTsr(tsrNew)
            endif
          enddo
        endif
        !
        if (lasttime) then
          kdone = 0
          do j=1,nts
            if (times(j) == badtime) kdone = kdone + 1
          enddo
          if (kdone == nts) exit loop1
        endif
        !
        ! Get another set of TS records
        lastCount = 0
        do j=1,nts
          chd => GetChdFromList(this%TvChdList, j)
          if (times(j) > smallest) then
            tsr => chd%timeSeries%GetCurrentTimeSeriesRecord()
          else
            tsr => chd%timeSeries%GetNextTimeSeriesRecord()
          endif
          if (associated(tsr)) then
            times(j) = tsr%tsrTime
            values(j) = tsr%tsrValue
            if (tsr%tsrTime == endtime) then
              lasttime = .true.
            endif
          else
            times(j) = badtime
            values(j) = DNODATA
            lastCount = lastCount + 1
          endif
        enddo
        call get_smallest(nts, times, ismall, smallest, allequal)
        !
      enddo loop1
    endif
    !
    deallocate(times)
    deallocate(values)
    return
  end subroutine InflateTimeSeries

  subroutine PopulateIchdvar(this, igrid)
    implicit none
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: i, iu, j, k, kper, n
    double precision :: ehead, shead
    !
    iu = this%IuOrig
    !
    ! Initialize ichdvar to 0 (not constant head)
    this%ichdvar = 0
    !
    ! Read through entire CHD input file (including all stress periods) to:
    !   Set ichdvar = -kper where CHD ehead = shead for all periods;
    !   Set ichdvar = kper where CHD ehead /= shead for any period.
    ! ABS(Ichdvar) indicates stress period in which CHD is first specified as active.
    !
    call GWF2CHD7AR(iu, igrid, .true.)   ! Store parameter data
    this%numtimeseries = 0
    perloop: do kper=1,NPER
      call GWF2CHD7RP(iu, igrid)
      chdloop: do n=1,NCHDS
        k = CHDS(1,n)  ! layer index (original model)
        i = CHDS(2,n)  ! row index
        j = CHDS(3,n)  ! column index
        shead = CHDS(4,n)
        ehead = CHDS(5,n)
        if (ehead /= shead) then
          if (this%ichdvar(j,i,k) <= 0) then
            if (this%ichdvar(j,i,k) < 0) then
              this%ichdvar(j,i,k) = -this%ichdvar(j,i,k)
            else
              this%ichdvar(j,i,k) = kper
            endif
            this%numtimeseries = this%numtimeseries + 1
          endif
        else
          if (this%ichdvar(j,i,k) == 0) then
            this%ichdvar(j,i,k) = -kper
          endif
        endif
      enddo chdloop
    enddo perloop
    ! Rewind CHD input file and deallocate CHD package variables
    rewind iu
    call GWF2CHD7DA(igrid)
    !
    ! Count all cells that are CHD in any stress period
    this%maxchds = 0
    do k=1,NLAY
      do i=1,NROW
        do j=1,NCOL
          if (this%ichdvar(j,i,k) /= 0) then
            this%maxchds = this%maxchds + 1
          endif
        enddo
      enddo
    enddo
    !
    return
  end subroutine PopulateIchdvar

  subroutine ReduceTimeSeries(this)
    ! dummy
    class(ChdPackageWriterType) :: this
    ! local
    integer :: i, j, nts
    type(ChdType), pointer :: chdi, chdj
    !
    nts = this%TvChdList%Count()
    if (nts > 1) then
      iloop: do i=2,nts
        chdi => GetChdFromList(this%TvChdList, i)
        jloop: do j=1,i-1
          chdj => GetChdFromList(this%TvChdList, j)
          if (associated(chdi%timeSeries) .and. associated(chdj%timeSeries)) then
            if (SameTimeSeries(chdi%timeSeries, chdj%timeSeries)) then
              chdi%AlternateTsName = chdj%timeSeries%Name
              call chdi%timeSeries%Clear(.true.)
              deallocate(chdi%timeSeries)
              cycle iloop
            endif
          endif
        enddo jloop
      enddo iloop
    endif
    !
    return
  end subroutine ReduceTimeSeries

  subroutine RewriteOutputFile(this)
    ! Rewrite output generated from CHD input to insert lines
    ! into Options block to list all required time-series files.
    ! dummy
    class(ChdPackageWriterType) :: this
    ! local
    integer :: i, istat, iu, j, linelen, maxb, nlines
    integer :: icol, istart, istop, idum
    double precision :: rdum
    logical :: inperblock
    character(len=MAXCHARLEN) :: blankline, line, word, msg
    type(LineListType), pointer :: lineList
    type(ChdType), pointer :: chd
    ! format
    10 format(a)
    20 format(2x,'TS6 FILEIN',2x,a)
    30 format(2x,a,2x,i0)
    !
    ! Copy contents of output file into lineList
    iu = this%fileobj%IUnit
    allocate(lineList)
    call lineList%InitializeLineList()
    rewind(iu)
    do
      read(iu,10,iostat=istat)line
      if (istat /= 0) exit
      call lineList%AddLine(line)
    enddo
    !
    ! Copy lines back into output file, with insertion of TIMESERIESFILE options
    blankline = ' '
    rewind(iu)
    maxb = this%GetMaxbound()
    nlines = lineList%CountLines()
    inperblock = .false.
    do i=1,nlines
      call lineList%GetLine(i,line)
      linelen = len_trim(line)
      write(iu,10)blankline(1:linelen)
      backspace(iu)
      if (line == 'END Options') then
        write(iu,20)trim(this%TsFile)
      endif
      icol = 1
      call URWORD(line, icol, istart, istop, 1, idum, rdum, 0, 0)
      word = line(istart:istop)
      select case (word)
      case ('MAXBOUND')
        write(iu,30)'MAXBOUND', maxb
      case ('AUXILIARY')
        msg = 'Auxiliary variables are not supported for conversion' // &
              ' of CHD Package.'
        call store_note(msg)
      case default
        if (line == 'END PERIOD') then
          inperblock = .false.
        endif
        if (inperblock) then
          call this%SubstituteTsName(line)
          write(iu,10)trim(line)
        else
          write(iu,10)trim(line)
          if (line(1:12) == 'BEGIN PERIOD') then
            inperblock = .true.
          endif
        endif
      end select
    enddo
    !
    call lineList%Clear(.true.)
    deallocate(lineList)
    !
    return
  end subroutine RewriteOutputFile

  subroutine SubstituteTsName(this, line)
    implicit none
    ! dummy
    class(ChdPackageWriterType), intent(inout) :: this
    character(len=*), intent(inout) :: line
    ! local
    integer :: i, ii, j, k, nchd
    character(len=LENTIMESERIESNAME) :: word
    type(ChdType), pointer :: chd
    !
    read(line,*)k,i,j,word
    if (word(1:4) /= 'CHTS') return
    !
    nchd = this%TvChdList%Count()
    chdloop: do ii=1,nchd
      chd => GetChdFromList(this%TvChdList, ii)
      if (associated(chd)) then
        if (chd%OriginalTsName == word) then
          if (chd%AlternateTsName /= ' ') then
            ! Rewrite line, substituting alternate ts name
            line = ' '
            write(line,this%fmat)k, i, j, trim(chd%AlternateTsName)
          endif
          exit chdloop
        endif
      endif
    enddo chdloop
    !
    return
  end subroutine SubstituteTsName

  subroutine WriteConstantListData(this, lineList)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    !
    integer :: i, ii, j, k, kold, kper1, n
    character(len=MAXCHARLEN) :: line, msg
    ! format
    20 format(2x,i0,2x,i0,2x,i0,2x,g15.8)
    10 format(a,' boundary removed at (',i0,',',i0,',',i0,')')
    !
    ! Write constant CHDs from CHD input
    if (this%IuOrig > 0) then
      do kold=1,NLAY
        k = this%Layptr(kold)
        do i=1,NROW
          do j=1,NCOL
            if (this%ichdvar(j,i,kold) < 0) then
              kper1 = -this%ichdvar(j,i,kold)
              if (this%kpercurrent >= kper1) then
                if (associated(masteridomain)) then
                  if (masteridomain(j,i,k)==0) then
                    if (verbose) then
                      msg = ''
                      write(msg,10)trim(this%PkgType),k,i,j
                      write(*,*)trim(msg)
                    endif
                    cycle
                  endif
                endif
                write(line,20)k,i,j,this%chdshead(j,i,kold)
                call lineList%AddLine(line)
              endif
            endif
          enddo
        enddo
      enddo
    endif
    !
    return
  end subroutine WriteConstantListData

  subroutine WriteIboundChdListData(this, lineList)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    integer :: i, ii, j, k, n, nibchd
    character(len=MAXCHARLEN) :: line, msg
    type(ChdType), pointer :: chd => null()
    ! format
    20 format(2x,i0,2x,i0,2x,i0,2x,g15.8)
    10 format(a,' boundary removed at (',i0,',',i0,',',i0,')')
    !
    ! Write CHDs from IBOUND
    n = this%IbChdList%Count()
    nibchd = n
    do ii=1,n
      chd => GetChdFromList(this%IbChdList, ii)
      j = chd%jcol
      i = chd%irow
      k = chd%klay
      if (associated(masteridomain)) then
        if (masteridomain(j,i,k)==0) then
          if (verbose) then
            msg = ''
            write(msg,10)trim(this%PkgType),k,i,j
            write(*,*)trim(msg)
          endif
          cycle
        endif
      endif
      write(line,20)k,i,j,chd%head
      call lineList%AddLine(line)
    enddo
    !
    return
  end subroutine WriteIboundChdListData

  subroutine WriteTimeSeriesFile(this)
    ! dummy
    class(ChdPackageWriterType) :: this
    ! local
    integer :: i, istrt, istp, iu, knodata, nts
    double precision :: time
    character(len=1), parameter :: dot='.', space=' '
    character(len=MAXCHARLEN) :: msg, tsfilename
    character(len=LENHUGELINE) :: allnames, tsrline
    character(len=LENTIMESERIESNAME) :: tsname
    integer, parameter      :: FIELDWID = 22
    character(len=FIELDWID) :: ctime, cvalue
    character(len=*), parameter :: fmtv = 'g22.15'
    character(len=*), parameter :: fmtpv = '(' // fmtv // ')'
    character(len=*), parameter :: fmtvxv = '(' // fmtv // ',1x,' // fmtv // ')'
    logical :: timewritten
    type(ChdType), pointer :: chd
    type(TimeSeriesRecordType), pointer :: tsr
    ! formats
    1 format()
    10 format(a,1x,a)
    20 format(2x,a,1x,a)
    !
    ! Ensure that all time series contain records for same set of times
    call this%InflateTimeSeries()
    ! Eliminate redundant time series
    call this%ReduceTimeSeries()
    !
    allnames = ' '
    !
    ! Open one file to hold all time series needed for CHD package
    tsfilename = trim(this%ModelBasename) // '_chd.ts'
    iu = GetUnit()
    call openfile(iu, 0, tsfilename, 'TS', filstat_opt='REPLACE')
    this%TsFile = tsfilename
    !
    nts = this%TvChdList%Count()
    !
    ! Write Attributes block
    ! Use METHOD LINEAREND to reproduce Leake & Prudic TWRI 6A2 behavior
    write(iu,10)'BEGIN','ATTRIBUTES'
    ! Build list of time-series names and write to NAMES attribute
    do i=1,nts
      chd => GetChdFromList(this%TvChdList,i)
      if (associated(chd%timeSeries)) then
        tsname = chd%timeSeries%Name
        if (allnames == ' ') then
          allnames = adjustl(tsname)
        else
          allnames = trim(allnames) // space // adjustl(tsname)
        endif
        call chd%timeSeries%Reset()
      endif
    enddo
    write(iu,20)'NAMES',trim(allnames)
    write(iu,20)'METHOD LINEAREND'
    write(iu,10)'END','ATTRIBUTES'
    !
    write(iu,1)
    write(iu,10)'BEGIN','TIMESERIES'
    ! Write all time-series records for all chd items
    timeloop: do
      ! Loop through CHDs, but skip if not time-varying or redundant
      knodata = 0
      tsrline = ' '
      timewritten = .false.
      chdloop: do i=1,nts
        chd => GetChdFromList(this%TvChdList,i)
        if (associated(chd%timeSeries)) then
          tsr => chd%timeSeries%GetNextTimeSeriesRecord()
          if (associated(tsr)) then
            if (.not. timewritten) then
              write(ctime,fmtpv)tsr%tsrTime
              tsrline = adjustl(ctime)
              write(cvalue,fmtpv)tsr%tsrValue
              tsrline = trim(tsrline) // space // adjustl(cvalue)
              !write(tsrline,fmtvxv)tsr%tsrTime, tsr%tsrValue
              time = tsr%tsrTime
              timewritten = .true.
              !istp = 2 * FIELDWID + 1
            else
              if (tsr%tsrTime /= time) then
                msg = 'Programmer error in WriteTimeSeriesFiles (CHD): time conflict.'
                call store_error(msg)
                call ustop()
              endif
              !write(iu,40,advance=adv)tsr%tsrValue
              write(cvalue,fmtpv)tsr%tsrValue
!              istrt = istp + 2
!              istp = istrt + FIELDWID - 1
!              tsrline(istrt:istp) = cvalue
              tsrline = trim(tsrline) // space // adjustl(cvalue)
            endif
          else
            knodata = knodata + 1
          endif
        else
          knodata = knodata + 1
        endif
        if (knodata == nts) exit timeloop
      enddo chdloop
      if (tsrline /= '') then
        write(iu,'(a)')trim(tsrline)
      endif
    enddo timeloop
    !
    write(iu,10)'END','TIMESERIES'
    !
    close(iu)
    !
    return
  end subroutine WriteTimeSeriesFile

  subroutine WriteTimeVaryingListData(this, lineList)
    implicit none
    ! dummy
    class(ChdPackageWriterType) :: this
    type(LineListType), pointer :: lineList
    ! local
    integer :: i, ii, j, k, n, nibchd
    character(len=MAXCHARLEN) :: line, msg
    type(ChdType), pointer :: chd => null()
    ! format
    20 format(2x,i0,2x,i0,2x,i0,2x,g15.8)
    10 format(a,' boundary removed at (',i0,',',i0,',',i0,')')
    !
    if (this%IuOrig > 0) then
      ! Write CHDs read from CHD package input file
      n = this%TvChdList%Count()
      do ii=1,n
        chd => GetChdFromList(this%TvChdList, ii)
        if (.not. chd%ActiveByStressPeriod(this%kpercurrent)) cycle
        j = chd%jcol
        i = chd%irow
        k = chd%klay
        if (associated(masteridomain)) then
          if (masteridomain(j,i,k)==0) then
            if (verbose) then
              msg = ''
              write(msg,10)trim(this%PkgType),k,i,j
              write(*,*)trim(msg)
            endif
            cycle
          endif
        else
          if (associated(this%Idomain)) then
            if (this%Idomain(j,i,k)==0) then
              if (verbose) then
                msg = ''
                write(msg,10)trim(this%PkgType),k,i,j
                write(*,*)trim(msg)
              endif
              cycle
            endif
          endif
        endif
        line = ' '
        if (chd%AlternateTsName == ' ') then
          write(line,this%fmat)k, i, j, trim(chd%timeSeries%Name)
        else
          write(line,this%fmat)k, i, j, trim(chd%AlternateTsName)
        endif
        call lineList%AddLine(line)
      enddo
    endif
    !
    return
  end subroutine WriteTimeVaryingListData

  ! Non-type-bound procedure

  subroutine get_smallest(idm, arr, ismall, smallest, allequal)
    ! dummy
    integer,                          intent(in)  :: idm
    double precision, dimension(idm), intent(in)  :: arr
    integer,                          intent(out) :: ismall
    double precision,                 intent(out) :: smallest
    logical,                          intent(out) :: allequal
    ! local
    integer :: i
    double precision ::val1
    !
    ismall = 0
    smallest = huge(smallest)
    allequal = .true.
    do i=1,idm
      if (arr(i) < smallest) then
        smallest = arr(i)
        ismall = i
      endif
      if (i==1) then
        val1 = arr(i)
      else
        if (arr(i) /= val1) allequal = .false.
      endif
    enddo
    !
    return
  end subroutine get_smallest

end module ChdPackageWriterModule
