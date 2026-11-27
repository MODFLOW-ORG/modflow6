module VsfSoilModelModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENMEMPATH, DZERO, DONE, DEP20
  use SimModule, only: store_error, count_errors, ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value
  use GwfVsfInputModule, only: GwfVsfParamFoundType

  implicit none
  private

  ! TODO_VSF: split this file in separate models

  public :: SoilModelType

  type, abstract :: SoilModelType
    integer(I4B) :: id !< the unique soil model id
    character(len=LENMEMPATH) :: input_mem_path !< the path for loading from input context
    character(len=LENMEMPATH) :: mem_path !< memory manager path for storing variables
    integer(I4B) :: nodes !< the number of nodes, dimensioning the model arrays
    integer(I4B), dimension(:), pointer, contiguous :: map => null() !< a map for node reduction
    real(DP), dimension(:), pointer, contiguous :: porosity => null() !< the model porosity per node
    real(DP), dimension(:), pointer, contiguous :: sat_res => null() !< the model residual saturation per node
  contains
    procedure :: create
    procedure(load_if), deferred :: load
    procedure :: destroy
    procedure :: destroy_base
    procedure(saturation_if), deferred :: saturation
    procedure(capacity_if), deferred :: capacity
    procedure(krelative_if), deferred :: krelative
    procedure :: kr_dkrdh
    procedure :: sat_capacity
    procedure :: pressure_head
  end type SoilModelType

  !> shift in pressure head for the numerical derivative of k_r
  real(DP), parameter, public :: DKR_SHIFT = 1.0e-6_DP

  abstract interface
    subroutine load_if(this, porosity, sat_res)
      import SoilModelType, DP
      class(SoilModelType), intent(inout) :: this
      real(DP), dimension(:), pointer, contiguous :: porosity
      real(DP), dimension(:), pointer, contiguous :: sat_res
    end subroutine
    function saturation_if(this, psi, i) result(s)
      import SoilModelType, DP, I4B
      class(SoilModelType), intent(inout) :: this
      real(DP) :: psi
      integer(I4B) :: i
      real(DP) :: s
    end function
    function capacity_if(this, psi, i) result(Cm)
      import SoilModelType, DP, I4B
      class(SoilModelType), intent(inout) :: this
      real(DP) :: psi
      integer(I4B) :: i
      real(DP) :: Cm
    end function
    function krelative_if(this, psi, i) result(kr)
      import SoilModelType, DP, I4B
      class(SoilModelType), intent(inout) :: this
      real(DP) :: psi
      integer(I4B) :: i
      real(DP) :: kr
    end function
  end interface

contains

  subroutine create(this, input_mem_path, mem_path, nodes, nodemap)
    class(SoilModelType), intent(inout) :: this
    character(len=*) :: input_mem_path
    character(len=*) :: mem_path
    integer(I4B) :: nodes
    integer(I4B), dimension(:), pointer, contiguous :: nodemap

    this%input_mem_path = input_mem_path
    this%mem_path = mem_path
    this%nodes = nodes
    this%map => nodemap

  end subroutine create

  subroutine destroy(this)
    class(SoilModelType), intent(inout) :: this

    call this%destroy_base()

  end subroutine destroy

  !> @brief Relative permeability and its derivative to pressure head
  !!
  !! Default: central difference of krelative. Soil models override this
  !! with a fused analytical evaluation.
  !<
  subroutine kr_dkrdh(this, psi, i, kr, dkr)
    class(SoilModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: kr !< relative permeability
    real(DP), intent(out) :: dkr !< d(kr)/d(psi)

    kr = this%krelative(psi, i)
    dkr = (this%krelative(psi + DKR_SHIFT, i) - &
           this%krelative(psi - DKR_SHIFT, i)) / (2.0_DP * DKR_SHIFT)

  end subroutine kr_dkrdh

  !> @brief Saturation and moisture capacity
  !!
  !! Default: separate calls. Soil models override this with a fused
  !! evaluation.
  !<
  subroutine sat_capacity(this, psi, i, s, cm)
    class(SoilModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: s !< saturation
    real(DP), intent(out) :: cm !< moisture capacity

    s = this%saturation(psi, i)
    cm = this%capacity(psi, i)

  end subroutine sat_capacity

  !> @brief Pressure head for a given effective saturation
  !!
  !! Inverse of the retention curve: returns psi such that the effective
  !! saturation (s - sr) / (1 - sr) at psi equals se. Returns 0 for se >= 1
  !! (soil models with an air-entry head return that head) and -1e20
  !! for se <= 0. Default: bisection on log10(-psi) between -1e-10 and
  !! -1e10. Soil models override this with the closed-form inverse.
  !<
  function pressure_head(this, se, i) result(psi)
    class(SoilModelType), intent(inout) :: this
    real(DP), intent(in) :: se !< effective saturation
    integer(I4B), intent(in) :: i !< node number
    real(DP) :: psi !< pressure head
    ! local
    real(DP) :: ulo, uhi, umid, sr, semid
    integer(I4B) :: it

    if (se >= DONE) then
      psi = DZERO
      return
    else if (se <= DZERO) then
      psi = -DEP20
      return
    end if
    sr = this%sat_res(i)
    ! effective saturation decreases with u = log10(-psi)
    ulo = -10.0_DP
    uhi = 10.0_DP
    do it = 1, 100
      umid = 0.5_DP * (ulo + uhi)
      semid = (this%saturation(-10.0_DP**umid, i) - sr) / (DONE - sr)
      if (semid > se) then
        ulo = umid
      else
        uhi = umid
      end if
    end do
    psi = -10.0_DP**(0.5_DP * (ulo + uhi))

  end function pressure_head

  subroutine destroy_base(this)
    class(SoilModelType), intent(inout) :: this

    this%map => null()
    this%porosity => null()
    this%sat_res => null()

  end subroutine destroy_base

end module VsfSoilModelModule
