module VsfHaverkampModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE
  use SimModule, only: store_error, count_errors, ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value
  use VsfSoilModelModule, only: SoilModelType, DKR_SHIFT
  use GwfVsfInputModule, only: GwfVsfParamFoundType
  implicit none
  private

  public :: HaverkampModelType

  ! TODO_VSF: document these parameters based on the tech ref
  type, extends(SoilModelType) :: HaverkampModelType
    real(DP), dimension(:), pointer, contiguous :: alpha => null()
    real(DP), dimension(:), pointer, contiguous :: n => null()
    real(DP), dimension(:), pointer, contiguous :: beta => null()
    real(DP), dimension(:), pointer, contiguous :: k => null()
  contains
    procedure :: load => load_haverkamp
    procedure :: saturation => saturation_haverkamp
    procedure :: capacity => capacity_haverkamp
    procedure :: krelative => krelative_haverkamp
    procedure :: kr_dkrdh => kr_dkrdh_haverkamp
    procedure :: sat_capacity => sat_capacity_haverkamp
    procedure :: destroy => destroy_haverkamp
  end type HaverkampModelType

contains

  !> @brief Load the Haverkamp model from input data
  !<
  subroutine load_haverkamp(this, porosity, sat_res)
    class(HaverkampModelType), intent(inout) :: this !< this instance
    real(DP), dimension(:), pointer, contiguous :: porosity !< porosity
    real(DP), dimension(:), pointer, contiguous :: sat_res !< residual saturation
    ! local
    type(GwfVsfParamFoundType) :: found

    this%porosity => porosity
    this%sat_res => sat_res

    call mem_allocate(this%alpha, this%nodes, 'ALPHAHVK', this%mem_path)
    call mem_allocate(this%n, this%nodes, 'NHVK', this%mem_path)
    call mem_allocate(this%beta, this%nodes, 'BETAHVK', this%mem_path)
    call mem_allocate(this%k, this%nodes, 'KHVK', this%mem_path)

    ! load alpha and n
    call mem_set_value(this%alpha, 'ALPHAHVK', this%input_mem_path, &
                       this%map, found%alphahvk)
    call mem_set_value(this%n, 'NHVK', this%input_mem_path, &
                       this%map, found%nhvk)
    call mem_set_value(this%beta, 'BETAHVK', this%input_mem_path, &
                       this%map, found%betahvk)
    call mem_set_value(this%k, 'KHVK', this%input_mem_path, &
                       this%map, found%khvk)
    if (.not. found%alphahvk) then
      call store_error('Haverkamp soil data (ALPHAHVK) &
        &not specified in input', .false.)
    end if
    if (.not. found%nhvk) then
      call store_error('Haverkamp soil data (NHVK) &
        &not specified in input', .false.)
    end if
    if (.not. found%betahvk) then
      call store_error('Haverkamp soil data (BETAHVK) &
        &not specified in input', .false.)
    end if
    if (.not. found%khvk) then
      call store_error('Haverkamp soil data (KHVK) &
        &not specified in input', .false.)
    end if
    if (count_errors() > 0) call ustop()

  end subroutine load_haverkamp

  function saturation_haverkamp(this, psi, i) result(s)
    class(HaverkampModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: s !< the saturation
    ! local
    real(DP) :: s_eff

    if (psi > DZERO) then
      s = DONE
    else
      s_eff = DONE / (DONE + (-this%alpha(i) * psi)**this%n(i))
      s = (1 - this%sat_res(i)) * s_eff + this%sat_res(i)
    end if

  end function saturation_haverkamp

  function capacity_haverkamp(this, psi, i) result(Cm)
    class(HaverkampModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: Cm !< the saturation
    ! local
    real(DP) :: num, denom

    if (psi > DZERO) then
      Cm = DZERO
    else
      num = this%n(i) * this%alpha(i) * &
            ((-this%alpha(i) * psi)**(this%n(i) - DONE))
      denom = (DONE + (-this%alpha(i) * psi)**this%n(i))**2
      Cm = this%porosity(i) * (DONE - this%sat_res(i)) * num / denom
    end if

  end function capacity_haverkamp

  function krelative_haverkamp(this, psi, i) result(kr)
    class(HaverkampModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: kr !< the saturation
    ! local

    if (psi > DZERO) then
      kr = DONE
    else
      kr = DONE / (DONE + (-this%beta(i) * psi)**this%k(i))
    end if

  end function krelative_haverkamp

  !> @brief Fused relative permeability and its analytical derivative
  !!
  !! Close to saturation the central difference of the default
  !! implementation is kept.
  !<
  subroutine kr_dkrdh_haverkamp(this, psi, i, kr, dkr)
    class(HaverkampModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: kr !< relative permeability
    real(DP), intent(out) :: dkr !< d(kr)/d(psi)
    ! local
    real(DP) :: x

    if (psi > -DKR_SHIFT) then
      kr = this%krelative(psi, i)
      dkr = (this%krelative(psi + DKR_SHIFT, i) - &
             this%krelative(psi - DKR_SHIFT, i)) / (2.0_DP * DKR_SHIFT)
      return
    end if

    x = (-this%beta(i) * psi)**this%k(i)
    kr = DONE / (DONE + x)
    dkr = this%k(i) * x / (-psi) * kr * kr

  end subroutine kr_dkrdh_haverkamp

  !> @brief Fused saturation and moisture capacity
  !<
  subroutine sat_capacity_haverkamp(this, psi, i, s, cm)
    class(HaverkampModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: s !< saturation
    real(DP), intent(out) :: cm !< moisture capacity
    ! local
    real(DP) :: x, se

    if (psi >= DZERO) then
      s = DONE
      cm = DZERO
      return
    end if

    x = (-this%alpha(i) * psi)**this%n(i)
    se = DONE / (DONE + x)
    s = (DONE - this%sat_res(i)) * se + this%sat_res(i)
    cm = this%porosity(i) * (DONE - this%sat_res(i)) * &
         this%n(i) * x / (-psi) * se * se

  end subroutine sat_capacity_haverkamp

  !> @brief clean up
  !<
  subroutine destroy_haverkamp(this)
    class(HaverkampModelType), intent(inout) :: this !< this instance

    this%porosity => null()
    this%sat_res => null()

    if (associated(this%alpha)) then
      call mem_deallocate(this%alpha)
      call mem_deallocate(this%n)
      call mem_deallocate(this%beta)
      call mem_deallocate(this%k)
    end if

    call this%destroy_base()

  end subroutine destroy_haverkamp

end module VsfHaverkampModule
