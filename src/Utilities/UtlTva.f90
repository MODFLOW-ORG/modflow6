!> @brief This module contains the UtlTva Module
!!
!! This module contains the code for accessing a generic set of time varying
!! grid arrays. The arrays are named by the user through the AUXILIARY option
!! and supplied per stress period through the Input Data Model, which also
!! handles any time-array series linked to them. A package that wants
!! time-varying gridded input declares utl-tva as a subpackage, creates a
!! UtlTvaType, and asks it for an array by name.
!<
module UtlTvaModule

  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: LENPACKAGENAME, LENMODELNAME, &
                             LENAUXNAME, LENMEMPATH, DZERO
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_deallocate, &
                                 mem_setptr, get_isize
  use MemoryManagerExtModule, only: mem_set_value
  use CharacterStringModule, only: CharacterStringType
  use BaseDisModule, only: DisBaseType

  implicit none
  private
  public :: UtlTvaType

  !> @brief Derived type for accessing time varying arrays
  !!
  !! Provides name-based access to the auxiliary arrays that the Input Data
  !! Model loads for this subpackage. The arrays themselves live in the input
  !! memory path and are refreshed there as the simulation advances through
  !! stress periods.
  !<
  type :: UtlTvaType

    character(len=LENMODELNAME) :: component = '' !< name of the parent component (typically model name)
    character(len=LENPACKAGENAME) :: subcomponent = '' !< name of the subcomponent (typically package name)
    character(len=LENMEMPATH) :: memoryPath = '' !< location in the memory manager for this object
    character(len=LENMEMPATH) :: input_mempath = '' !< location of the input loaded by the Input Data Model
    character(len=LENAUXNAME), dimension(:), pointer, &
      contiguous :: auxname => null() !< vector of auxiliary variable names
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< auxiliary names as loaded from input
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< the loaded arrays, dimensioned (nvar, isize)
    integer(I4B), pointer :: iout => null() !< unit number for output
    integer(I4B), pointer :: nvar => null() !< number of variables
    integer(I4B), pointer :: isize => null() !< length of each array
    integer(I4B), pointer :: iprpak => null() !< flag for printing input
    class(DisBaseType), pointer :: dis => null() !< model discretization object

  contains

    procedure :: initialize
    procedure :: has
    procedure :: fill
    procedure :: tva_da
    procedure, private :: allocate_scalars
    procedure, private :: source_options

  end type UtlTvaType

contains

  !> @ brief Initialize the TVA type
  !!
  !! Set up the object against the input that the Input Data Model has loaded
  !< for this subpackage, and point at the loaded arrays.
  subroutine initialize(this, dis, isize, input_mempath, iout, component, &
                        subcomponent)
    ! dummy
    class(UtlTvaType) :: this !<  UtlTvaType
    class(DisBaseType), pointer, intent(in) :: dis !<  discretization package
    integer(I4B), intent(in) :: isize !<  size of array
    character(len=*), intent(in) :: input_mempath !<  memory path of the loaded input
    integer(I4B), intent(in) :: iout !<  unit number for output
    character(len=*), intent(in) :: component !<  character string containing component
    character(len=*), intent(in) :: subcomponent !<  character string containing subcomponent

    ! construct the memory paths
    this%component = component
    this%subcomponent = subcomponent
    this%memoryPath = create_mem_path(this%component, &
                                      trim(this%subcomponent)//'-TVA')
    this%input_mempath = input_mempath

    ! allocate scalar variables
    call this%allocate_scalars()

    ! assign member values
    this%iout = iout
    this%isize = isize

    ! set pointers
    this%dis => dis

    ! source options from the loaded input
    call this%source_options()
  end subroutine initialize

  !> @ Check to see if variable was included by user
  function has(this, auxname) result(found)
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object
    character(len=*), intent(in) :: auxname
    ! return
    logical(LGP) :: found
    ! local
    integer(I4B) :: i

    found = .false.
    do i = 1, this%nvar
      if (this%auxname(i) == auxname) then
        found = .true.
        exit
      end if
    end do
  end function has

  !> @ brief Allocate scalars
  !<
  subroutine allocate_scalars(this)
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object

    ! allocate scalars in memory manager
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%nvar, 'NVAR', this%memoryPath)
    call mem_allocate(this%isize, 'ISIZE', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%auxname, LENAUXNAME, 0, 'AUXNAME', this%memoryPath)
    call mem_allocate(this%auxname_cst, LENAUXNAME, 0, 'AUXNAME_CST', &
                      this%memoryPath)

    ! initialize
    this%iout = 0
    this%nvar = 0
    this%isize = 0
    this%iprpak = 0
  end subroutine allocate_scalars

  !> @ brief Source options from the loaded input
  !!
  !! Read the auxiliary names and the print flag from the input memory path,
  !! then point at the array of loaded values. Any time-array series attached
  !< to these arrays is handled by the Input Data Model.
  subroutine source_options(this)
    ! dummy
    class(UtlTvaType) :: this
    ! local
    integer(I4B) :: n
    integer(I4B) :: isize_aux
    logical(LGP) :: found
    ! formats
    character(len=*), parameter :: fmtiprpak = &
      &"(4x,'TVA INFORMATION WILL BE PRINTED TO LISTING FILE.')"

    ! print flag
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, found)
    if (found .and. this%iprpak /= 0) write (this%iout, fmtiprpak)

    ! auxiliary variable names
    call get_isize('AUXILIARY', this%input_mempath, isize_aux)
    if (isize_aux > 0) then
      this%nvar = isize_aux
      call mem_reallocate(this%auxname, LENAUXNAME, this%nvar, &
                          'AUXNAME', this%memoryPath)
      call mem_reallocate(this%auxname_cst, LENAUXNAME, this%nvar, &
                          'AUXNAME_CST', this%memoryPath)
      call mem_set_value(this%auxname_cst, 'AUXILIARY', this%input_mempath, &
                         found, release=.false.)
      do n = 1, this%nvar
        this%auxname(n) = this%auxname_cst(n)
      end do
    end if

    ! point at the loaded arrays; they are not present until a period block
    ! supplies values for at least one variable
    call get_isize('AUXVAR', this%input_mempath, isize_aux)
    if (isize_aux > 0) then
      call mem_setptr(this%auxvar, 'AUXVAR', this%input_mempath)
    end if
  end subroutine source_options

  !> @ brief Fill the passed array with values for the named variable
  !<
  subroutine fill(this, dblvec, auxname)
    ! dummy
    class(UtlTvaType) :: this !< UtlTvaType object
    character(len=*), intent(in) :: auxname
    real(DP), dimension(:), intent(inout) :: dblvec
    ! local
    integer(I4B) :: i
    integer(I4B) :: ivar
    integer(I4B) :: nu
    integer(I4B) :: nr

    ivar = 0
    do i = 1, this%nvar
      if (this%auxname(i) == auxname) then
        ivar = i
        exit
      end if
    end do

    if (ivar == 0) then
      write (errmsg, '(a)') 'Variable not found in TVA input: '//trim(auxname)
      call store_error(errmsg, terminate=.TRUE.)
    end if

    ! nothing has been supplied for any variable yet
    if (.not. associated(this%auxvar)) return

    ! copy values from the loaded array into the passed array
    do nr = 1, this%dis%nodes
      nu = this%dis%get_nodeuser(nr)
      dblvec(nr) = this%auxvar(ivar, nu)
    end do
  end subroutine fill

  !> @ brief Deallocate
  !<
  subroutine tva_da(this)
    ! dummy variables
    class(UtlTvaType) :: this !< UtlTvaType object

    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxname_cst, 'AUXNAME_CST', this%memoryPath)
    call mem_deallocate(this%iout)
    call mem_deallocate(this%nvar)
    call mem_deallocate(this%isize)
    call mem_deallocate(this%iprpak)

    nullify (this%auxvar)
    nullify (this%dis)
  end subroutine tva_da

end module UtlTvaModule
