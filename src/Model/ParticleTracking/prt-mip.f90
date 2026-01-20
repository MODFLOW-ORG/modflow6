module PrtMipModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LINELENGTH
  use NumericalPackageModule, only: NumericalPackageType
  use BlockParserModule, only: BlockParserType
  use BaseDisModule, only: DisBaseType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value, memorystore_remove
  use SimVariablesModule, only: idm_context
  use SimModule, only: store_error
  use PrtMipInputModule, only: PrtMipParamFoundType

  implicit none
  private
  public :: PrtMipType
  public :: mip_cr

  type, extends(NumericalPackageType) :: PrtMipType
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< aquifer porosity
    real(DP), dimension(:), pointer, contiguous :: retfactor => null() !< retardation factor
    integer(I4B), dimension(:), pointer, contiguous :: izone => null() !< zone number
  contains
    procedure :: mip_ar
    procedure :: mip_da
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
  end type PrtMipType

contains

  !> @brief Create a model input object
  subroutine mip_cr(mip, name_model, input_mempath, inunit, iout, dis)
    ! dummy
    type(PrtMipType), pointer :: mip
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    ! formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'MIP MODEL INPUT PACKAGE', &
       &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! Create the object
    allocate (mip)
    !
    ! Create name and memory path
    call mip%set_names(1, name_model, 'MIP', 'MIP', input_mempath)
    !
    ! Allocate scalars
    call mip%allocate_scalars()
    !
    ! Set variables
    mip%inunit = inunit
    mip%iout = iout
    !
    ! Set pointers
    mip%dis => dis
    !
    ! Print a message identifying the package if enabled
    if (inunit > 0) &
      write (iout, fmtheader) input_mempath

  end subroutine mip_cr

  !> @brief Deallocate memory
  subroutine mip_da(this)
    class(PrtMipType) :: this
    !
    ! Deallocate input memory
    call memorystore_remove(this%name_model, 'MIP', idm_context)
    !
    ! Deallocate parent package
    call this%NumericalPackageType%da()
    !
    ! Deallocate arrays
    call mem_deallocate(this%porosity)
    call mem_deallocate(this%retfactor)
    call mem_deallocate(this%izone)

  end subroutine mip_da

  subroutine allocate_scalars(this)
    class(PrtMipType) :: this
    call this%NumericalPackageType%allocate_scalars()
  end subroutine allocate_scalars

  !> @brief Allocate arrays
  subroutine allocate_arrays(this, nodes)
    class(PrtMipType) :: this
    integer(I4B), intent(in) :: nodes
    ! local
    integer(I4B) :: i
    !
    ! Allocate
    call mem_allocate(this%porosity, nodes, 'POROSITY', this%memoryPath)
    call mem_allocate(this%retfactor, nodes, 'RETFACTOR', this%memoryPath)
    call mem_allocate(this%izone, nodes, 'IZONE', this%memoryPath)
    !
    do i = 1, nodes
      this%porosity(i) = DZERO
      this%retfactor(i) = DONE
      this%izone(i) = 0
    end do

  end subroutine allocate_arrays

  !> @ brief Initialize package inputs
  subroutine mip_ar(this)
    ! dummy variables
    class(PrtMipType), intent(inout) :: this !< PrtMipType object
    ! local variables
    character(len=LINELENGTH) :: errmsg
    type(PrtMipParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map => null()
    !
    ! set map to convert user input data into reduced data
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser
    !
    ! Allocate arrays
    call this%allocate_arrays(this%dis%nodes)
    !
    ! Source array inputs from IDM
    call mem_set_value(this%porosity, 'POROSITY', this%input_mempath, &
                       map, found%porosity)
    call mem_set_value(this%retfactor, 'RETFACTOR', this%input_mempath, &
                       map, found%retfactor)
    call mem_set_value(this%izone, 'IZONE', this%input_mempath, map, &
                       found%izone)
    !
    ! Ensure POROSITY was found
    if (.not. found%porosity) then
      write (errmsg, '(a)') 'Error in GRIDDATA block: POROSITY not found'
      call store_error(errmsg)
    end if

  end subroutine mip_ar

end module PrtMipModule
