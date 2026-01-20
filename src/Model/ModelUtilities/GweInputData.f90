module GweInputDataModule

  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, LENMEMPATH, DEP3

  implicit none
  private
  public :: GweInputDataType
  public :: gweshared_dat_cr
  public :: set_gwe_dat_ptrs

  !> Data for sharing among multiple packages.  Originally read in from
  !< the EST package

  type GweInputDataType

    ! dim
    integer(I4B) :: nnodes !< number of cells

    ! strings
    character(len=LENMEMPATH) :: memoryPath = '' !< the location in the memory manager where the variables are stored

    ! est data to be share across multiple packages
    real(DP), pointer :: gwerhow => null() !< Density of water (for GWE purposes, a constant scalar)
    real(DP), pointer :: gwecpw => null() !< Heat capacity of water (non-spatially varying)
    real(DP), pointer :: gwelatheatvap => null() !< latent heat of vaporization
    real(DP), dimension(:), pointer, contiguous :: gwerhos => null() !< Density of the aquifer material
    real(DP), dimension(:), pointer, contiguous :: gwecps => null() !< Heat capacity of solids (spatially varying)

  contains

    ! public
    procedure, public :: set_gwe_dat_ptrs
    procedure, public :: gweshared_dat_da
    ! private
    procedure, private :: set_gwe_scalar_ptrs
    procedure, private :: set_gwe_array_ptrs

  end type GweInputDataType

contains

  !> @brief Allocate the shared data
  !<
  subroutine gweshared_dat_cr(gwe_dat)
    ! modules
    ! dummy
    type(GweInputDataType), pointer :: gwe_dat !< the input data block

    ! create the object
    allocate (gwe_dat)

  end subroutine gweshared_dat_cr

  !> @brief Allocate and read data from EST
  !!
  !! EST data, including heat capacity of water (cpw), density of water
  !! (rhow), latent heat of vaporization (latheatvap), heat capacity of
  !! the aquifer material (cps), and density of the aquifer material
  !! (rhow) is used among other packages and is therefore stored in a
  !! separate class
  !<
  subroutine set_gwe_dat_ptrs(this, gwerhow, gwecpw, gwelatheatvap, &
                              gwerhos, gwecps)
    ! dummy
    class(GweInputDataType) :: this !< the input data block
    real(DP), intent(in), pointer :: gwerhow !< ptr to density of water specified in EST
    real(DP), intent(in), pointer :: gwecpw !< ptr to heat capacity of water specified in EST
    real(DP), intent(in), pointer :: gwelatheatvap !< ptr to latent heat of vaporization specified in EST
    real(DP), dimension(:), pointer, contiguous :: gwerhos !< ptr to sptially-variably density of aquifer material specified in EST
    real(DP), dimension(:), pointer, contiguous :: gwecps !< ptr to sptially-variably heat capacity of aquifer material specified in EST

    ! allocate scalars
    call this%set_gwe_scalar_ptrs(gwerhow, gwecpw, gwelatheatvap)

    ! allocate arrays
    call this%set_gwe_array_ptrs(gwerhos, gwecps)

  end subroutine set_gwe_dat_ptrs

  !> @brief Set pointers to scalars read by the EST package
  !! for use by other packages
  !!
  !! Set pointers to GWE-related scalars and arrays for use
  !! by multiple packages.  For example, a package capable of
  !! simulating evaporation will need access to latent heat of
  !! of vaporization.
  !<
  subroutine set_gwe_scalar_ptrs(this, gwerhow, gwecpw, gwelatheatvap)
    ! dummy
    class(GweInputDataType) :: this !< GweInputDataType object
    real(DP), pointer, intent(in) :: gwerhow !< density of water
    real(DP), pointer, intent(in) :: gwecpw !< mass-based heat capacity of water
    real(DP), pointer, intent(in), optional :: gwelatheatvap !< latent heat of vaporization

    ! set the pointers
    ! fixed density of water to be used by GWE (vs GWT-based density)
    this%gwerhow => gwerhow
    ! spatially constant heat capacity of water   ! kluge note: "specific heat" (which is heat capacity per unit mass) is probably the more correct term
    this%gwecpw => gwecpw
    ! latent heat of vaporization
    if (present(gwelatheatvap)) then
      this%gwelatheatvap => gwelatheatvap
    end if

  end subroutine set_gwe_scalar_ptrs

  !> @brief Set pointers to data arrays read by the EST package
  !! for use by other packages
  !!
  !! Set pointers to GWE-related arrays for use by multiple packages
  !<
  subroutine set_gwe_array_ptrs(this, gwerhos, gwecps)
    ! dummy
    class(GweInputDataType) :: this !< GweInputDataType object
    real(DP), dimension(:), pointer, contiguous, intent(in) :: gwerhos
    real(DP), dimension(:), pointer, contiguous, intent(in) :: gwecps

    ! set the pointers
    ! spatially-variable density of aquifer solid material
    this%gwerhos => gwerhos
    ! spatially-variable heat capacity of aquifer solid material
    this%gwecps => gwecps

  end subroutine set_gwe_array_ptrs

  !> @ brief Deallocate memory
  !!
  !!  Set pointers to null
  !<
  subroutine gweshared_dat_da(this)
    ! dummy
    class(GweInputDataType) :: this !< the input data block

    ! scalars
    this%gwelatheatvap => null()
    this%gwerhow => null()
    this%gwecpw => null()

    ! arrays
    this%gwerhos => null()
    this%gwecps => null()

  end subroutine gweshared_dat_da

end module GweInputDataModule
