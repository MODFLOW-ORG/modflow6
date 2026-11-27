!> @brief Brooks-Corey soil characteristic model
!!
!! With x = -alpha psi, alpha = 1 / |air-entry head| and the pore-size
!! distribution index lambda (Brooks and Corey, 1964):
!!
!!   se = 1                    for x <= 1 (saturated above the air entry)
!!   se = x**(-lambda)         for x > 1
!!   kr = se**(3 + 2/lambda)   (Burdine)
!!
!! The moisture capacity and the derivative of kr jump at the air-entry
!! head, a property of the model.
!<
module VsfBrooksCoreyModule
  use KindModule, only: I4B, DP
  use ConstantsModule, only: DZERO, DONE, DEP20
  use SimModule, only: store_error, count_errors, ustop
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryManagerExtModule, only: mem_set_value
  use VsfSoilModelModule, only: SoilModelType
  use GwfVsfInputModule, only: GwfVsfParamFoundType
  implicit none
  private

  public :: BrooksCoreyModelType

  type, extends(SoilModelType) :: BrooksCoreyModelType
    real(DP), dimension(:), pointer, contiguous :: alpha => null() !< inverse of the absolute air-entry head
    real(DP), dimension(:), pointer, contiguous :: lambda => null() !< pore-size distribution index
  contains
    procedure :: load => load_brookscorey
    procedure :: saturation => saturation_brookscorey
    procedure :: capacity => capacity_brookscorey
    procedure :: krelative => krelative_brookscorey
    procedure :: kr_dkrdh => kr_dkrdh_brookscorey
    procedure :: sat_capacity => sat_capacity_brookscorey
    procedure :: pressure_head => pressure_head_brookscorey
    procedure :: destroy => destroy_brookscorey
  end type BrooksCoreyModelType

contains

  !> @brief Load the Brooks-Corey model from input data
  !<
  subroutine load_brookscorey(this, porosity, sat_res)
    class(BrooksCoreyModelType), intent(inout) :: this !< this instance
    real(DP), dimension(:), pointer, contiguous :: porosity !< porosity
    real(DP), dimension(:), pointer, contiguous :: sat_res !< residual saturation
    ! local
    type(GwfVsfParamFoundType) :: found

    this%porosity => porosity
    this%sat_res => sat_res

    call mem_allocate(this%alpha, this%nodes, 'ALPHABC', this%mem_path)
    call mem_allocate(this%lambda, this%nodes, 'LAMBDABC', this%mem_path)

    call mem_set_value(this%alpha, 'ALPHABC', this%input_mem_path, &
                       this%map, found%alphabc)
    call mem_set_value(this%lambda, 'LAMBDABC', this%input_mem_path, &
                       this%map, found%lambdabc)
    if (.not. found%alphabc) then
      call store_error('Brooks-Corey soil data (ALPHABC) &
        &not specified in input', .false.)
    else if (any(this%alpha <= DZERO)) then
      call store_error('Brooks-Corey soil data (ALPHABC) &
        &must be positive', .false.)
    end if
    if (.not. found%lambdabc) then
      call store_error('Brooks-Corey soil data (LAMBDABC) &
        &not specified in input', .false.)
    else if (any(this%lambda <= DZERO)) then
      call store_error('Brooks-Corey soil data (LAMBDABC) &
        &must be positive', .false.)
    end if
    if (count_errors() > 0) call ustop()

  end subroutine load_brookscorey

  function saturation_brookscorey(this, psi, i) result(s)
    class(BrooksCoreyModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: s !< the saturation
    ! local
    real(DP) :: cm

    call this%sat_capacity(psi, i, s, cm)

  end function saturation_brookscorey

  function capacity_brookscorey(this, psi, i) result(cm)
    class(BrooksCoreyModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: cm !< the moisture capacity
    ! local
    real(DP) :: s

    call this%sat_capacity(psi, i, s, cm)

  end function capacity_brookscorey

  function krelative_brookscorey(this, psi, i) result(kr)
    class(BrooksCoreyModelType), intent(inout) :: this !< this instance
    real(DP) :: psi !< the pressure head in cell
    integer(I4B) :: i !< the cell number
    real(DP) :: kr !< the relative permeability
    ! local
    real(DP) :: dkr

    call this%kr_dkrdh(psi, i, kr, dkr)

  end function krelative_brookscorey

  !> @brief Relative permeability (Burdine) and its derivative
  !!
  !! kr = se**(3 + 2/lambda) = x**(-(2 + 3 lambda)), x = -alpha psi
  !<
  subroutine kr_dkrdh_brookscorey(this, psi, i, kr, dkr)
    class(BrooksCoreyModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: kr !< relative permeability
    real(DP), intent(out) :: dkr !< d(kr)/d(psi)
    ! local
    real(DP) :: x, eta

    x = -this%alpha(i) * psi
    if (x <= DONE) then
      kr = DONE
      dkr = DZERO
      return
    end if
    eta = 2.0_DP + 3.0_DP * this%lambda(i)
    kr = exp(-eta * log(x))
    dkr = eta * kr / (-psi)

  end subroutine kr_dkrdh_brookscorey

  !> @brief Saturation and moisture capacity
  !<
  subroutine sat_capacity_brookscorey(this, psi, i, s, cm)
    class(BrooksCoreyModelType), intent(inout) :: this
    real(DP), intent(in) :: psi !< pressure head
    integer(I4B), intent(in) :: i !< node number
    real(DP), intent(out) :: s !< saturation
    real(DP), intent(out) :: cm !< moisture capacity
    ! local
    real(DP) :: x, se

    x = -this%alpha(i) * psi
    if (x <= DONE) then
      s = DONE
      cm = DZERO
      return
    end if
    se = exp(-this%lambda(i) * log(x))
    s = (DONE - this%sat_res(i)) * se + this%sat_res(i)
    cm = this%porosity(i) * (DONE - this%sat_res(i)) * &
         this%lambda(i) * se / (-psi)

  end subroutine sat_capacity_brookscorey

  !> @brief Pressure head for a given effective saturation
  !!
  !! Inverse of se = (-alpha psi)**(-lambda). For se >= 1 the air-entry
  !! head -1/alpha is returned.
  !<
  function pressure_head_brookscorey(this, se, i) result(psi)
    class(BrooksCoreyModelType), intent(inout) :: this
    real(DP), intent(in) :: se !< effective saturation
    integer(I4B), intent(in) :: i !< node number
    real(DP) :: psi !< pressure head

    if (se >= DONE) then
      psi = -DONE / this%alpha(i)
    else if (se <= DZERO) then
      psi = -DEP20
    else
      psi = -exp(-log(se) / this%lambda(i)) / this%alpha(i)
    end if

  end function pressure_head_brookscorey

  !> @brief clean up
  !<
  subroutine destroy_brookscorey(this)
    class(BrooksCoreyModelType), intent(inout) :: this !< this instance

    this%porosity => null()
    this%sat_res => null()

    if (associated(this%alpha)) then
      call mem_deallocate(this%alpha)
      call mem_deallocate(this%lambda)
    end if

    call this%destroy_base()

  end subroutine destroy_brookscorey

end module VsfBrooksCoreyModule
