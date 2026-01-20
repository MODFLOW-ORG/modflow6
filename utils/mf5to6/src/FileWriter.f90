module FileWriterModule

  use ConstantsModule, only: LINELENGTH, LENFTYPE, MAXCHARLEN
  use FileTypeModule, only: FileType
  use GLOBAL, only: iout
  use InputOutputModule, only: GetUnit, openfile
  use SimPHMFModule, only: store_warning
  use UtilitiesModule, only: close_file

  private
  public :: FileWriterType, FileType

  type :: FileWriterType
    logical,                 public :: Active = .false.
    type(FileType), pointer, public :: fileobj => null()
    character(len=4),        public :: PkgType = ''
  contains
    procedure, public :: InitializeFile => initialize_by_name
    procedure, public :: CloseFile
  end type FileWriterType

contains

  subroutine initialize_by_name(this, fname, ftype, pkgname)
    implicit none
    ! dummy
    class(FileWriterType), intent(inout) :: this
    character(len=*), intent(in)    :: fname
    character(len=*), intent(in)    :: ftype
    character(len=*), intent(in), optional :: pkgname
    ! local
    character(len=7) :: fstatus = 'REPLACE'
    character(len=10) :: cdate, ctime
    character(len=2) :: dtrim, mtrim, htrim
    ! formats
    10 format('# ',a,' input file, prepared by MF5to6 on ', &
              a,' at ',a,'.')
    !
    ! Allocate the FileType member
    if (.not. associated(this%fileobj)) then
      allocate(this%fileobj)
      call this%fileobj%Initialize()
    endif
    !
    ! Store file info and open file
    this%fileobj%FName = fname
    this%fileobj%FType = ftype
    if (present(pkgname)) then
      this%fileobj%PkgName = pkgname
    endif
    !
    ! Open file only if file writer is active
    if (this%Active) then
      this%fileobj%IUnit = GetUnit()
      call openfile(this%fileobj%IUnit, iout, fname, ftype, 'FORMATTED', &
                    'SEQUENTIAL', fstatus)
      !
      ! Time stamp
      call date_and_time(cdate, ctime)
      ! trim month
      if (cdate(5:5)=='0') then
        mtrim = cdate(6:6)
      else
        mtrim = cdate(5:6)
      endif
      ! trim day
      if (cdate(7:7)=='0') then
        dtrim = cdate(8:8)
      else
        dtrim = cdate(7:8)
      endif
      cdate = trim(mtrim) // '/' // trim(dtrim) // '/' // cdate(1:4)
      ! trim hour
      if (ctime(1:1)=='0') then
        htrim = ctime(2:2)
      else
        htrim = ctime(1:2)
      endif
      ctime = trim(htrim) // ':' // ctime(3:4) // ':' // ctime(5:6)
      !
      write(this%fileobj%IUnit,10) trim(ftype), trim(cdate), trim(ctime)
    endif
    !
    return
  end subroutine initialize_by_name

  subroutine CloseFile(this)
    implicit none
    ! dummy
    class(FileWriterType) :: this
    ! local
    character(len=MAXCHARLEN) :: ermsg
    !
    if (close_file(this%fileobj%IUnit,100,0.1d0)) then
      this%fileobj%IUnit = 0
    else
      ermsg = 'Unable to close file: ' // trim(this%fileobj%FName)
      call store_warning(ermsg)
    endif
    !
    return
  end subroutine CloseFile

end module FileWriterModule
