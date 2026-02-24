module BaseModelModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMODELNAME, LINELENGTH, LENMEMPATH
  use ListModule, only: ListType
  use BaseDisModule, only: DisBaseType
  use InputOutputModule, only: openfile, getunit
  use VersionModule, only: write_listfile_header
  implicit none

  private
  public :: BaseModelType, CastAsBaseModelClass, AddBaseModelToList, &
            GetBaseModelFromList

  !> @brief Highest level model type. All models extend this parent type.
  type :: BaseModelType
    character(len=LENMEMPATH) :: memoryPath !< the location in the memory manager where the variables are stored
    character(len=LENMODELNAME), pointer :: name => null() !< name of the model
    character(len=LINELENGTH), pointer :: filename => null() !< input file name
    character(len=3), pointer :: macronym => null() !< 3 letter model acronym (GWF, GWT, ...)
    integer(I4B), pointer :: idsoln => null() !< id of the solution model is in
    integer(I4B), pointer :: id => null() !< model id
    integer(I4B), pointer :: iout => null() !< output unit number
    integer(I4B), pointer :: inewton => null() !< newton-raphson flag
    integer(I4B), pointer :: iprpak => null() !< integer flag to echo input
    integer(I4B), pointer :: iprflow => null() !< flag to print simulated flows
    integer(I4B), pointer :: ipakcb => null() !< save_flows flag
    integer(I4B), pointer :: nja => null() !< number of connections
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< ibound
    real(DP), dimension(:), pointer, contiguous :: flowja => null() !< intercell flows
    type(ListType), pointer :: bndlist => null() !< array of boundary packages
    class(DisBaseType), pointer :: dis => null() !< discretization object
  contains
    procedure :: model_df
    procedure :: model_ar
    procedure :: model_rp
    procedure :: model_dt
    procedure :: model_ot
    procedure :: model_fp
    procedure :: model_da
    procedure :: allocate_scalars
    procedure :: model_message
    procedure :: create_lstfile
  end type BaseModelType

contains

  !> @brief Define the model
  subroutine model_df(this)
    class(BaseModelType) :: this
  end subroutine model_df

  !> @brief Allocate and read
  subroutine model_ar(this)
    class(BaseModelType) :: this
  end subroutine model_ar

  !> @brief Read and prepare
  subroutine model_rp(this)
    class(BaseModelType) :: this
  end subroutine model_rp

  !> @brief Calculate time step length
  subroutine model_dt(this)
    class(BaseModelType) :: this
  end subroutine model_dt

  !> @brief Output results
  subroutine model_ot(this)
    class(BaseModelType) :: this
  end subroutine model_ot

  !> @brief Write line to model iout
  subroutine model_message(this, line, fmt)
    ! dummy
    class(BaseModelType) :: this
    character(len=*), intent(in) :: line
    character(len=*), intent(in), optional :: fmt
    ! local
    character(len=LINELENGTH) :: cfmt

    if (present(fmt)) then
      cfmt = fmt
    else
      cfmt = '(1x,a)'
    end if

    write (this%iout, trim(cfmt)) trim(line)
  end subroutine model_message

  !> @brief Final processing
  subroutine model_fp(this)
    class(BaseModelType) :: this
  end subroutine model_fp

  !> @brief Allocate scalar variables
  subroutine allocate_scalars(this, modelname)
    ! modules
    use MemoryManagerModule, only: mem_allocate
    ! dummy
    class(BaseModelType) :: this
    character(len=*), intent(in) :: modelname

    call mem_allocate(this%name, LENMODELNAME, 'NAME', this%memoryPath)
    call mem_allocate(this%macronym, 3, 'MACRONYM', this%memoryPath)
    call mem_allocate(this%id, 'ID', this%memoryPath)
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%inewton, 'INEWTON', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%iprflow, 'IPRFLOW', this%memoryPath)
    call mem_allocate(this%ipakcb, 'IPAKCB', this%memoryPath)
    call mem_allocate(this%idsoln, 'IDSOLN', this%memoryPath)

    this%name = modelname
    this%macronym = ''
    this%idsoln = 0
    this%id = 0
    this%iout = 0
    this%iprpak = 0
    this%iprflow = 0
    this%ipakcb = 0
    this%inewton = 0 !default is standard formulation
  end subroutine allocate_scalars

  !> @brief Deallocate
  subroutine model_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    ! dummy
    class(BaseModelType) :: this

    ! deallocate strings
    call mem_deallocate(this%name, 'NAME', this%memoryPath)
    call mem_deallocate(this%macronym, 'MACRONYM', this%memoryPath)
    !
    ! deallocate scalars
    call mem_deallocate(this%id)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%inewton)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%iprflow)
    call mem_deallocate(this%ipakcb)
    call mem_deallocate(this%idsoln)
  end subroutine model_da

  function CastAsBaseModelClass(obj) result(res)
    class(*), pointer, intent(inout) :: obj
    class(BaseModelType), pointer :: res

    res => null()
    if (.not. associated(obj)) return

    select type (obj)
    class is (BaseModelType)
      res => obj
    end select
  end function CastAsBaseModelClass

  subroutine AddBaseModelToList(list, model)
    ! dummy
    type(ListType), intent(inout) :: list
    class(BaseModelType), pointer, intent(inout) :: model
    ! local
    class(*), pointer :: obj

    obj => model
    call list%Add(obj)
  end subroutine AddBaseModelToList

  function GetBaseModelFromList(list, idx) result(res)
    ! dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(BaseModelType), pointer :: res
    ! local
    class(*), pointer :: obj

    obj => list%GetItem(idx)
    res => CastAsBaseModelClass(obj)
  end function GetBaseModelFromList

  subroutine create_lstfile(this, lst_fname, model_fname, defined, headertxt)
    ! dummy
    class(BaseModelType) :: this
    character(len=*), intent(inout) :: lst_fname
    character(len=*), intent(in) :: model_fname
    logical(LGP), intent(in) :: defined
    character(len=*), intent(in) :: headertxt
    ! local
    integer(I4B) :: i, istart, istop

    ! set list file name if not provided
    if (.not. defined) then

      ! initialize
      lst_fname = ' '
      istart = 0
      istop = len_trim(model_fname)

      ! identify '.' character position from back of string
      do i = istop, 1, -1
        if (model_fname(i:i) == '.') then
          istart = i
          exit
        end if
      end do

      ! if not found start from string end
      if (istart == 0) istart = istop + 1

      ! set list file name
      lst_fname = model_fname(1:istart)
      istop = istart + 3
      lst_fname(istart:istop) = '.lst'
    end if

    ! create the list file
    this%iout = getunit()
    call openfile(this%iout, 0, lst_fname, 'LIST', filstat_opt='REPLACE')

    ! write list file header
    call write_listfile_header(this%iout, headertxt)
  end subroutine create_lstfile

end module BaseModelModule
