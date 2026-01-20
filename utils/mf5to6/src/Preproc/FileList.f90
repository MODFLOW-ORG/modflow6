module FileListModule

  use ConstantsModule,     only: LINELENGTH, MAXCHARLEN
  use ConstantsPHMFModule, only: FCUNKNOWN, FCINPUT, FCDATAIN, FCDATABIN, &
                                 FCDATAOUT, FCOUTPUT, FCDATABOUT
  use FileTypeModule,      only: FileType, ConstructFileType, CastAsFileType
  use InputOutputModule,   only: same_word
  use ListModule,          only: ListType
  use ListNodeModule,      only: ListNodeType
  use SimPHMFModule,       only: store_warning, store_error, ustop
  use UtilitiesModule,     only: close_file

  private
  public :: FileListType

  type :: FileListType
    type(ListType), pointer, public :: files => null()
  contains
    ! Public procedures
    procedure, public  :: Initialize => initialize_list
    generic, public    :: AddFile => add_file, add_file_object
    generic, public    :: GetFile => get_file_by_index, get_file_by_type
    procedure, public  :: GetFileByUnit => get_file_by_unit
    procedure, public  :: GetNextFile
    procedure, public  :: CloseAll => close_all
    procedure, public  :: NCount
    procedure, public  :: ContainsUnit
    ! Private procedures
    procedure, private :: add_file
    procedure, private :: add_file_object
    procedure, private :: get_file_by_type
    procedure, private :: get_file_by_index
  end type FileListType

contains

  ! Methods of FileListType
  subroutine initialize_list(this)
    implicit none
    ! dummy arguments
    class(FileListType), intent(inout) :: this
    !
    if (associated(this%files)) then
      call this%files%Clear(.true.)
      deallocate(this%files)
    endif
    allocate(this%files)
    return
  end subroutine initialize_list

  subroutine add_file(this, fname, ftype, iu, FCode, pkgname)
    implicit none
    ! dummy arguments
    class(FileListType), intent(inout) :: this
    character(len=LINELENGTH), intent(in) :: fname
    character(len=*),    intent(in) :: ftype
    integer,             intent(in) :: iu
    integer,             intent(in) :: FCode
    character(len=*), optional, intent(in) :: pkgname
    ! local variables
    type(FileType), pointer :: newFile => null()
    character(len=100) :: msg
    !
    if (fname/='' .and. ftype/='' .and. iu>=0) then
      call ConstructFileType(newFile)
      if (FCode == FCDATAOUT .or. FCode == FCOUTPUT .or. &
          FCode == FCDATABOUT) then
        ! Don't allow slashes or dots in output file name.
        newFile%fname = localize_file_name(fname)
      else
        newFile%fname = fname
      endif
      newFile%ftype = ftype
      newFile%iunit = iu
      if (present(pkgname)) then
        newfile%PkgName = pkgname
      endif
      if (FCode/=FCUNKNOWN) then
        newFile%FCode = FCode
        if ((FCode==FCINPUT .or. FCode==FCDATAIN .or. FCode==FCDATABIN) &
            .and. iu<1) then
          msg = 'Error in FileListType%add_file: IU < 1 for an INPUT file' // &
                ' named "' // trim(fname) // '" of type: ' // trim(ftype)
          call store_error(msg)
          call ustop()
        endif
      else
        msg = 'Error in FileListType%add_file: FCode is UNKNOWN for file "' // &
              trim(fname) // '" of type: ' // trim(ftype)
        call store_error(msg)
        call ustop()
      endif
      call this%add_file_object(newFile)
    else
      msg = 'Error in FileListType%add_file: name, type, or unit is invalid'
      call store_error(msg)
      call ustop()
    endif
    !
    return
  end subroutine add_file

  subroutine add_file_object(this, fileobj)
    ! dummy
    class(FileListType), intent(inout) :: this
    type(FileType), pointer, intent(inout) :: fileobj
    class(*), pointer :: obj => null()
    !
    obj => fileobj
    call this%files%Add(obj)
    !
    return
  end subroutine add_file_object

  logical function ContainsUnit(this, iu)
    ! dummy
    class(FileListType), intent(inout) :: this
    integer, intent(in) :: iu
    ! local
    integer :: i, n
    type(FileType), pointer :: fil => null()
    !
    ContainsUnit = .false.
    n = this%files%Count()
    do i=1,n
      fil => this%get_file_by_index(i)
      if (fil%IUnit == iu) then
        ContainsUnit = .true.
        exit
      endif
    enddo
    !
    return
  end function

  function get_file_by_unit(this, iu, killonfailure) result(file)
    ! If killonfailure is absent and file is not found,
    ! error is written and processing stops.
    implicit none
    ! dummy arguments
    class(FileListType), intent(inout) :: this
    integer, intent(in) :: iu
    logical, optional, intent(in) :: killonfailure
    type(FileType), pointer :: file
    ! local variables
    character(len=500) :: msg
    class(*), pointer :: obj => null()
    type(FileType), pointer :: filePtr => null()
    logical :: killonfailurelocal
    ! format
    10 format('Error: No file found associated with file unit: ',i0)
    !
    if (present(killonfailure)) then
      killonfailurelocal = killonfailure
    else
      killonfailurelocal = .true.
    endif
    !
    file => null()
    if (associated(this%files)) then
      call this%files%Reset()
      obj => this%files%GetNextItem()
      do while(associated(obj))
        filePtr => CastAsFileType(obj)
        if (filePtr%IUnit==iu) then
          file => filePtr
          exit
        endif
        obj => this%files%GetNextItem()
      enddo
    endif
    !
    if (.not. associated(file)) then
      if (killonfailurelocal) then
        write(msg,10)iu
        call store_error(msg)
        call ustop()
      endif
    endif
    !
    return
  end function get_file_by_unit

  function get_file_by_type(this, ftype) result(file)
    implicit none
    ! dummy arguments
    class(FileListType), intent(inout) :: this
    character(len=*),    intent(in) :: ftype
    type(FileType), pointer :: file
    ! local variables
    class(*), pointer :: obj => null()
    type(FileType), pointer :: filePtr => null()
    character(len=500) :: msg
    ! format
    10 format('Error: No file found of type: ',a)
    !
    file => null()
    call this%files%Reset()
    obj => this%files%GetNextItem()
    do while(associated(obj))
      filePtr => CastAsFileType(obj)
      if (same_word(ftype, filePtr%ftype)) then
        file => filePtr
        exit
      endif
      obj => this%files%GetNextItem()
    enddo
    !
    if (.not. associated(file)) then
      write(msg,10)trim(ftype)
      call store_error(ftype)
      call ustop()
    endif
    !
    return
  end function get_file_by_type

  function get_file_by_index(this, indx) result(fil)
    implicit none
    ! dummy arguments
    class(FileListType), intent(inout) :: this
    integer, intent(in) :: indx
    type(FileType), pointer :: fil
    ! local variables
    class(*), pointer :: obj => null()
    !
    fil => null()
    obj => this%files%GetItem(indx)
    if (associated(obj)) then
      fil => CastAsFileType(obj)
    endif
    return
  end function get_file_by_index

  function GetNextFile(this) result(fil)
    implicit none
    ! dummy
    class(FileListType), intent(inout) :: this
    type(FileType), pointer :: fil
    ! local
    class(*), pointer :: obj => null()
    !
    fil => null()
    obj => this%files%GetNextItem()
    if (associated(obj)) then
      fil => CastAsFileType(obj)
    endif
    !
    return
  end function GetNextFile

  subroutine close_all(this)
    implicit none
    ! dummy arguments
    class(FileListType), intent(inout) :: this
    ! local variables
    type(FileType), pointer :: fileptr
    integer :: iu
    class(*), pointer :: obj => null()
    character(len=MAXCHARLEN) :: ermsg, ermsgio
    !
    obj => this%files%GetItem(1)
    fileptr => CastAsFileType(obj)
    do while (associated(fileptr))
      iu = fileptr%iunit
      if (iu > 0) then
        if (close_file(iu,100,0.1d0)) then
          fileptr%IUnit = 0
        else
          ermsg = 'Unable to close file: ' // trim(fileptr%FName)
          call store_warning(ermsg)
          call store_warning(ermsgio)
        endif
      endif
      obj => this%files%GetNextItem()
      fileptr => CastAsFileType(obj)
    enddo
    !
    return
  end subroutine close_all

  integer function NCount(this)
    implicit none
    ! dummy
    class(FileListType), intent(inout) :: this
    !
    NCount = this%files%Count()
    !
    return
  end function NCount

  function localize_file_name(fnamein) result(fnameout)
    ! Eliminate dots; change / or \ to _
    implicit none
    ! dummy
    character(len=*), intent(in) :: fnamein
    character(len=LINELENGTH)    :: fnameout
    ! local
    integer :: i, j, leng
    !logical :: relative
    !
    fnameout = ' '
    leng = len_trim(fnamein)
    j = 0
    do i=1,leng
      if (fnamein(i:i) == '.') then
        if (i<leng) then
          if (fnamein(i+1:i+1) == '.' .or. fnamein(i+1:i+1) == '/' &
              .or. fnamein(i+1:i+1) == '\') then
            ! Dot character indicates a relative path
            cycle
          endif
        endif
      endif
      j = j + 1
      if (fnamein(i:i) == '\' .or. fnamein(i:i) == '/') then
        fnameout(j:j) = '_'
      else
        fnameout(j:j) = fnamein(i:i)
      endif
    enddo
    !
    return
  end function localize_file_name

end module FileListModule

