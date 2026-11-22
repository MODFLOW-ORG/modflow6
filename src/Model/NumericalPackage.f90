!> @brief This module contains the base numerical package type
!!
!! This module contains the base model package class that is extended
!! by all model packages.
!!
!<
module NumericalPackageModule
  ! -- modules
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, &
                             LENMEMPATH, LENFTYPE, LINELENGTH, &
                             LENVARNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  use MemoryHelperModule, only: create_mem_path

  implicit none
  private
  public NumericalPackageType

  type :: NumericalPackageType

    ! -- strings
    character(len=LENMODELNAME) :: name_model = '' !< the name of the model that contains this package
    character(len=LENPACKAGENAME) :: packName = '' !< name of the package
    character(len=LENMEMPATH) :: memoryPath = '' !< the location in the memory manager where the variables are stored
    character(len=LENMEMPATH) :: memoryPathModel = '' !< the location in the memory manager where the variables
                                                                                 !! of the parent model are stored
    character(len=LENMEMPATH) :: input_mempath = '' !< input context mempath
    character(len=LINELENGTH), pointer :: input_fname => null() !< input file name
    character(len=LENFTYPE) :: filtyp = '' !< file type (CHD, DRN, RIV, etc.)
    character(len=LENFTYPE), pointer :: package_type => null() !< package type (same as filtyp) stored in memory manager

    ! -- integers
    integer(I4B), pointer :: id => null() !< consecutive package number in model
    integer(I4B), pointer :: inunit => null() !< unit number for input file
    integer(I4B), pointer :: iout => null() !< unit number for writing package output
    integer(I4B), pointer :: inewton => null() !< newton flag
    integer(I4B), pointer :: iasym => null() !< package causes matrix asymmetry
    integer(I4B), pointer :: iprpak => null() !< integer flag to echo input
    integer(I4B), pointer :: iprflow => null() !< flag to print simulated flows
    integer(I4B), pointer :: ipakcb => null() !< output flows (-1, 0, 1) - save_flows
    integer(I4B), pointer :: ionper => null() !< stress period for next data
    integer(I4B), pointer :: lastonper => null() !< last value of ionper (for checking)
    !
    ! -- derived types
    type(BlockParserType) :: parser !< parser object for reading blocks of information
    class(DisBaseType), pointer :: dis => null() !< model discretization object

  contains
    procedure :: set_names
    procedure :: allocate_scalars
    procedure :: da
    procedure :: read_check_ionper
    procedure :: get_block_data
  end type NumericalPackageType
  !
contains
  !
  !> @brief Set package names
  !!
  !!  Method to assign the filtyp (ftype), the model name, and package name for
  !!  a package. This method also creates the memoryPath and memoryPathModel that
  !!  is used by the memory manager when variables are allocated.
  !!
  !<
  subroutine set_names(this, ibcnum, name_model, pakname, ftype, input_mempath)
    ! -- dummy variables
    class(NumericalPackageType), intent(inout) :: this !< NumericalPackageType object
    integer(I4B), intent(in) :: ibcnum !< unique package number
    character(len=*), intent(in) :: name_model !< name of the model
    character(len=*), intent(in) :: pakname !< name of the package
    character(len=*), intent(in) :: ftype !< package type
    character(len=*), optional, intent(in) :: input_mempath !< input_mempath
    !
    ! -- set names
    this%filtyp = ftype
    this%name_model = name_model
    if (present(input_mempath)) this%input_mempath = input_mempath
    if (pakname == '') then
      write (this%packName, '(a, i0)') trim(ftype)//'-', ibcnum
    else
      !
      ! -- Ensure pakname has no spaces
      if (index(trim(pakname), ' ') > 0) then
        errmsg = 'Package name contains spaces: '//trim(pakname)
        call store_error(errmsg)
        errmsg = 'Remove spaces from name.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
      !
      this%packName = pakname
    end if
    this%memoryPath = create_mem_path(name_model, this%packName)
    this%memoryPathModel = create_mem_path(name_model)
  end subroutine set_names

  !> @brief Allocate package scalars
  !!
  !!  Allocate and initialize base numerical package scalars.
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy variables
    class(NumericalPackageType) :: this !< NumericalPackageType object
    ! -- local variables
    integer(I4B), pointer :: imodelnewton => null()
    integer(I4B), pointer :: imodelprpak => null()
    integer(I4B), pointer :: imodelprflow => null()
    integer(I4B), pointer :: imodelpakcb => null()
    logical(LGP) :: found
    !
    ! -- allocate scalars
    call mem_allocate(this%input_fname, LINELENGTH, 'INPUT_FNAME', &
                      this%memoryPath)
    call mem_allocate(this%package_type, LENFTYPE, 'PACKAGE_TYPE', &
                      this%memoryPath)
    call mem_allocate(this%id, 'ID', this%memoryPath)
    call mem_allocate(this%inunit, 'INUNIT', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%inewton, 'INEWTON', this%memoryPath)
    call mem_allocate(this%iasym, 'IASYM', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%iprflow, 'IPRFLOW', this%memoryPath)
    call mem_allocate(this%ipakcb, 'IPAKCB', this%memoryPath)
    !
    call mem_allocate(this%ionper, 'IONPER', this%memoryPath)
    call mem_allocate(this%lastonper, 'LASTONPER', this%memoryPath)
    !
    ! -- set pointer to model variables
    call mem_setptr(imodelnewton, 'INEWTON', this%memoryPathModel)
    call mem_setptr(imodelprpak, 'IPRPAK', this%memoryPathModel)
    call mem_setptr(imodelprflow, 'IPRFLOW', this%memoryPathModel)
    call mem_setptr(imodelpakcb, 'IPAKCB', this%memoryPathModel)
    !
    ! -- initialize
    this%input_fname = ''
    this%package_type = this%filtyp
    this%id = 0
    this%inunit = 0
    this%iout = 0
    this%inewton = imodelnewton
    this%iasym = 0
    this%iprpak = imodelprpak
    this%iprflow = imodelprflow
    this%ipakcb = imodelpakcb
    this%ionper = 0
    this%lastonper = 0
    !
    ! -- nullify unneeded pointers
    imodelnewton => null()
    imodelprpak => null()
    imodelprflow => null()
    imodelpakcb => null()
    !
    ! -- update input filename
    if (this%input_mempath /= '') then
      call mem_set_value(this%input_fname, 'INPUT_FNAME', &
                         this%input_mempath, found)
    end if
  end subroutine allocate_scalars

  !> @brief Deallocate package scalars
  !!
  !!  Deallocate and initialize base numerical package scalars.
  !!
  !<
  subroutine da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(NumericalPackageType) :: this !< NumericalPackageType object
    !
    ! -- deallocate
    call mem_deallocate(this%input_fname, 'INPUT_FNAME', this%memoryPath)
    call mem_deallocate(this%package_type, 'PACKAGE_TYPE', this%memoryPath)
    call mem_deallocate(this%id)
    call mem_deallocate(this%inunit)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%inewton)
    call mem_deallocate(this%iasym)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%iprflow)
    call mem_deallocate(this%ipakcb)
    call mem_deallocate(this%ionper)
    call mem_deallocate(this%lastonper)
  end subroutine da

  !> @brief Check ionper
  !!
  !!  Generic method to read and check ionperiod, which is used to determine
  !!  if new period data should be read from the input file. The check of
  !!  ionperiod also makes sure periods are increasing in subsequent period
  !!  data blocks.
  !!
  !<
  subroutine read_check_ionper(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy variables
    class(NumericalPackageType), intent(inout) :: this !< NumericalPackageType object
    !
    ! -- save last value and read period number
    this%lastonper = this%ionper
    this%ionper = this%parser%GetInteger()
    !
    ! -- make check
    if (this%ionper <= this%lastonper) then
      write (errmsg, '(a, i0, a, i0, a, i0, a)') &
        'Error in stress period ', kper, &
        '. Period numbers not increasing.  Found ', this%ionper, &
        ' but last period block was assigned ', this%lastonper, '.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
  end subroutine read_check_ionper

  !> @brief Read griddata block for a package
  !!
  !!  Generic method to read data in the GRIDDATA block for a package.
  !!
  !<
  subroutine get_block_data(this, tags, lfound, varinames)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy variables
    class(NumericalPackageType) :: this !< NumericalPackageType object
    character(len=24), dimension(:), intent(in) :: tags !< vector with variable tags
    logical, dimension(:), intent(inout) :: lfound !< boolean vector indicating of a variable tag was found
    character(len=24), dimension(:), intent(in), optional :: varinames !< optional vector of variable names
    ! -- local variables
    logical :: lkeyword
    logical :: endOfBlock
    integer(I4B) :: nsize
    integer(I4B) :: j
    character(len=24) :: tmpvar
    character(len=LENVARNAME) :: varname
    character(len=LINELENGTH) :: keyword
    character(len=:), allocatable :: line
    integer(I4B) :: istart, istop, lloc
    integer(I4B), dimension(:), pointer, contiguous :: aint
    real(DP), dimension(:), pointer, contiguous :: adbl
    !
    ! -- initialize nsize
    nsize = size(tags)
    do
      call this%parser%GetNextLine(endOfBlock)
      if (endOfBlock) exit
      call this%parser%GetStringCaps(keyword)
      call this%parser%GetRemainingLine(line)
      lkeyword = .false.
      lloc = 1
      tag_iter: do j = 1, nsize
        if (trim(adjustl(keyword)) == trim(adjustl(tags(j)))) then
          lkeyword = .true.
          lfound(j) = .true.
          if (present(varinames)) then
            tmpvar = adjustl(varinames(j))
          else
            tmpvar = adjustl(tags(j))
          end if
          varname = tmpvar(1:LENVARNAME)
          if (keyword(1:1) == 'I') then
            call mem_setptr(aint, trim(varname), trim(this%memoryPath))
            call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                          this%parser%iuactive, aint, tags(j))
          else
            call mem_setptr(adbl, trim(varname), trim(this%memoryPath))
            call this%dis%read_grid_array(line, lloc, istart, istop, this%iout, &
                                          this%parser%iuactive, adbl, tags(j))
          end if
          exit tag_iter
        end if
      end do tag_iter
      if (.not. lkeyword) then
        write (errmsg, '(a,a)') 'Unknown GRIDDATA tag: ', &
          trim(keyword)
        call store_error(errmsg)
        call this%parser%StoreErrorUnit()
      end if
    end do
  end subroutine get_block_data

end module NumericalPackageModule

