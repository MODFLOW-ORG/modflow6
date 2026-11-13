!> @brief Time-Varying Grid (TVG) Utility
!!
!! Reads one or more named full-grid arrays that can change by stress period.
!! A parent package declares a TVG6 file in its options block, names the
!! auxiliary variables it expects, and calls fill() to copy the current
!! period's values into an array of its own.
!!
!! Input is loaded by the Input Data Model grid-array loader, which stores
!! values compressed: only cells whose value is not DNODATA are kept, with
!! NODEULIST giving the user node number of each stored entry. fill()
!! expands that back onto the caller's reduced-node vector.
!<
module UtlTvgModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LENAUXNAME, LENPACKAGENAME, &
                             LENCOMPONENTNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use MemoryManagerModule, only: mem_allocate, mem_reallocate, mem_deallocate, &
                                 mem_setptr, get_isize
  use MemoryManagerExtModule, only: mem_set_value
  use MemoryHelperModule, only: create_mem_path
  use CharacterStringModule, only: CharacterStringType
  use BaseDisModule, only: DisBaseType
  use TdisModule, only: kper

  implicit none
  private
  public :: UtlTvgType

  type :: UtlTvgType
    character(len=LENCOMPONENTNAME) :: component = '' !< name of the component
    character(len=LENPACKAGENAME) :: subcomponent = '' !< name of the subcomponent (typically package name)
    character(len=LENMEMPATH) :: memoryPath = '' !< location in the memory manager for this object
    character(len=LENMEMPATH) :: input_mempath = '' !< location of the input loaded by the Input Data Model
    character(len=LENAUXNAME), dimension(:), pointer, &
      contiguous :: auxname => null() !< vector of auxiliary variable names
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< auxiliary names as loaded from input
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< loaded values, dimensioned (nvar, nbound)
    integer(I4B), dimension(:), pointer, contiguous :: nodeulist => null() !< user node number of each stored entry
    integer(I4B), pointer :: nbound => null() !< number of stored entries in the current period
    class(DisBaseType), pointer :: dis => null() !< discretization package
    integer(I4B), pointer :: iout => null() !< unit number for output
    integer(I4B), pointer :: nvar => null() !< number of variables
    integer(I4B), pointer :: iprpak => null() !< print input flag
  contains
    procedure :: initialize
    procedure :: has
    procedure :: fill
    procedure :: tvg_da
    procedure, private :: allocate_scalars
    procedure, private :: source_options
  end type UtlTvgType

contains

  !> @brief Initialize the TVG type
  !<
  subroutine initialize(this, dis, input_mempath, iout, component, subcomponent)
    ! dummy
    class(UtlTvgType) :: this !< UtlTvgType object
    class(DisBaseType), pointer, intent(in) :: dis !< discretization package
    character(len=*), intent(in) :: input_mempath !< memory path of the loaded input
    integer(I4B), intent(in) :: iout !< unit number for output
    character(len=*), intent(in) :: component !< character string containing component
    character(len=*), intent(in) :: subcomponent !< character string containing subcomponent

    ! construct the memory paths
    this%component = component
    this%subcomponent = subcomponent
    this%memoryPath = create_mem_path(this%component, &
                                      trim(this%subcomponent)//'-TVG')
    this%input_mempath = input_mempath

    call this%allocate_scalars()

    this%iout = iout
    this%dis => dis

    call this%source_options()
  end subroutine initialize

  !> @brief Check to see if variable was included by user
  !<
  function has(this, auxname) result(found)
    ! dummy
    class(UtlTvgType) :: this !< UtlTvgType object
    character(len=*), intent(in) :: auxname !< name to look for
    ! return
    logical(LGP) :: found
    ! local
    integer(I4B) :: n

    found = .false.
    do n = 1, this%nvar
      if (this%auxname(n) == auxname) then
        found = .true.
        exit
      end if
    end do
  end function has

  !> @brief Allocate scalars
  !<
  subroutine allocate_scalars(this)
    ! dummy
    class(UtlTvgType) :: this !< UtlTvgType object

    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%nvar, 'NVAR', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%auxname, LENAUXNAME, 0, 'AUXNAME', this%memoryPath)
    call mem_allocate(this%auxname_cst, LENAUXNAME, 0, 'AUXNAME_CST', &
                      this%memoryPath)

    this%iout = 0
    this%nvar = 0
    this%iprpak = 0
  end subroutine allocate_scalars

  !> @brief Source options from the loaded input
  !<
  subroutine source_options(this)
    ! dummy
    class(UtlTvgType) :: this !< UtlTvgType object
    ! local
    integer(I4B) :: n
    integer(I4B) :: isize
    logical(LGP) :: found
    ! formats
    character(len=*), parameter :: fmtiprpak = &
      &"(4x,'TVG INFORMATION WILL BE PRINTED TO LISTING FILE.')"

    ! print flag; the definition file renames PRINT_INPUT to IPRPAK, so the
    ! loaded variable is found under the internal name
    call mem_set_value(this%iprpak, 'IPRPAK', this%input_mempath, found)
    if (found .and. this%iprpak /= 0) write (this%iout, fmtiprpak)

    ! auxiliary variable names
    call get_isize('AUXILIARY', this%input_mempath, isize)
    if (isize > 0) then
      this%nvar = isize
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

  end subroutine source_options

  !> @brief Fill the passed array with values for the named variable
  !!
  !! The grid-array loader stores only the cells it read, so values are
  !! expanded through NODEULIST. Cells the loader did not supply are left
  !! as the caller had them, which lets a parent package keep its own
  !! defaults or the values from an earlier period.
  !<
  subroutine fill(this, dblvec, auxname)
    ! dummy
    class(UtlTvgType) :: this !< UtlTvgType object
    real(DP), dimension(:), intent(inout) :: dblvec !< array to fill
    character(len=*), intent(in) :: auxname !< name of the variable to copy
    ! local
    integer(I4B) :: i
    integer(I4B) :: ivar
    integer(I4B) :: nu
    integer(I4B) :: nr
    integer(I4B) :: isize
    integer(I4B), pointer :: iper

    ivar = 0
    do i = 1, this%nvar
      if (this%auxname(i) == auxname) then
        ivar = i
        exit
      end if
    end do

    if (ivar == 0) then
      write (errmsg, '(a)') 'Variable not found in TVG input: '//trim(auxname)
      call store_error(errmsg, terminate=.TRUE.)
    end if

    ! Resolve the loader's arrays on every call rather than caching them at
    ! initialize: they do not exist until a period block has been read, and
    ! the loader may reallocate them as later periods are loaded.
    call get_isize('AUXVAR', this%input_mempath, isize)
    if (isize <= 0) return
    call mem_setptr(this%auxvar, 'AUXVAR', this%input_mempath)
    call mem_setptr(this%nodeulist, 'NODEULIST', this%input_mempath)
    call mem_setptr(this%nbound, 'NBOUND', this%input_mempath)

    ! Only apply values the loader read for the current stress period. When a
    ! period block is absent the loader leaves the previous period's data in
    ! place and IPER no longer matches, so the caller keeps what it had.
    call mem_setptr(iper, 'IPER', this%input_mempath)
    if (iper /= kper) return

    ! expand the stored entries onto the caller's reduced-node vector
    do i = 1, this%nbound
      nu = this%nodeulist(i)
      nr = this%dis%get_nodenumber(nu, 0)
      if (nr > 0) then
        dblvec(nr) = this%auxvar(ivar, i)
      end if
    end do
  end subroutine fill

  !> @brief Deallocate
  !<
  subroutine tvg_da(this)
    ! dummy
    class(UtlTvgType) :: this !< UtlTvgType object

    call mem_deallocate(this%iout)
    call mem_deallocate(this%nvar)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxname_cst, 'AUXNAME_CST', this%memoryPath)

    nullify (this%auxvar)
    nullify (this%nodeulist)
    nullify (this%nbound)
    nullify (this%dis)
  end subroutine tvg_da

end module UtlTvgModule
