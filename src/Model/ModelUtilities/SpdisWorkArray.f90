module SpdisWorkArrayModule
  use ConstantsModule, only: DZERO
  use KindModule, only: I4B, DP, LGP
  implicit none
  private

  !> Helper class with work arrays for the SPDIS calculation in NPF
  !<
  type, public :: SpdisWorkArrayType
    real(DP), allocatable, dimension(:) :: vi
    real(DP), allocatable, dimension(:) :: di
    real(DP), allocatable, dimension(:) :: viz
    real(DP), allocatable, dimension(:) :: diz
    real(DP), allocatable, dimension(:) :: nix
    real(DP), allocatable, dimension(:) :: niy
    real(DP), allocatable, dimension(:) :: wix
    real(DP), allocatable, dimension(:) :: wiy
    real(DP), allocatable, dimension(:) :: wiz
    real(DP), allocatable, dimension(:) :: bix
    real(DP), allocatable, dimension(:) :: biy
  contains
    procedure :: create
    procedure :: is_created
    procedure :: reset
    procedure :: destroy
  end type SpdisWorkArrayType

contains

  !< @brief Create the worker arrays
  !<
  subroutine create(this, nr_conns)
    class(SpdisWorkArrayType) :: this
    integer(I4B) :: nr_conns

    allocate (this%vi(nr_conns))
    allocate (this%di(nr_conns))
    allocate (this%viz(nr_conns))
    allocate (this%diz(nr_conns))
    allocate (this%nix(nr_conns))
    allocate (this%niy(nr_conns))
    allocate (this%wix(nr_conns))
    allocate (this%wiy(nr_conns))
    allocate (this%wiz(nr_conns))
    allocate (this%bix(nr_conns))
    allocate (this%biy(nr_conns))

  end subroutine create

  !< @brief True when created (= allocated)
  !<
  function is_created(this) result(created)
    class(SpdisWorkArrayType) :: this
    logical(LGP) :: created

    created = allocated(this%vi)

  end function is_created

  !< @brief Sets all arrays to zero
  !<
  subroutine reset(this)
    class(SpdisWorkArrayType) :: this

    this%vi(:) = DZERO
    this%di(:) = DZERO
    this%viz(:) = DZERO
    this%diz(:) = DZERO
    this%nix(:) = DZERO
    this%niy(:) = DZERO
    this%wix(:) = DZERO
    this%wiy(:) = DZERO
    this%wiz(:) = DZERO
    this%bix(:) = DZERO
    this%biy(:) = DZERO

  end subroutine reset

  !< @brief Clean up memory (only when create has been called)
  !<
  subroutine destroy(this)
    class(SpdisWorkArrayType) :: this

    deallocate (this%vi)
    deallocate (this%di)
    deallocate (this%viz)
    deallocate (this%diz)
    deallocate (this%nix)
    deallocate (this%niy)
    deallocate (this%wix)
    deallocate (this%wiy)
    deallocate (this%wiz)
    deallocate (this%bix)
    deallocate (this%biy)

  end subroutine destroy

end module SpdisWorkArrayModule
