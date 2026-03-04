module TspExchangeMoverModule
  use TspMvtModule

  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENVARNAME
  use VirtualModelModule, only: VirtualModelType
  use TspFmiModule, only: TspFmiType
  implicit none
  private

  public :: xmvt_cr

  type, public, extends(TspMvtType) :: TspExchangeMoverType
    class(VirtualModelType), pointer :: model1 => null() !< virtual model 1
    class(VirtualModelType), pointer :: model2 => null() !< virtual model 2
    real(DP), dimension(:), pointer, contiguous :: flow => null() !< the MVR flow
    real(DP), dimension(:), pointer, contiguous :: qty => null() !< concentration, temp, ...
  contains
    procedure :: mvt_rp => xmvt_rp
    procedure :: xmvt_cf
    procedure :: mvt_fc => xmvt_fc    
  end type TspExchangeMoverType

contains
  
  subroutine xmvt_cr(mvt, name_model, inunit, iout, fmi1, eqnsclfac, &
                    depvartype, gwfmodelname1, gwfmodelname2, fmi2)
    type(TspExchangeMoverType), pointer :: mvt
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi1
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    character(len=LENVARNAME), intent(in) :: depvartype !< dependent variable type ('concentration' or 'temperature')
    character(len=*), intent(in), optional :: gwfmodelname1
    character(len=*), intent(in), optional :: gwfmodelname2
    type(TspFmiType), intent(in), target, optional :: fmi2

    allocate (mvt)
    call mvt%mvt_init(name_model, inunit, iout, fmi1, &
                      eqnsclfac, depvartype, gwfmodelname1, &
                      gwfmodelname2, fmi2)

  end subroutine xmvt_cr

  !> @brief Read and prepare mover transport object  
  !<
  subroutine xmvt_rp(this)
    class(TspExchangeMoverType) :: this
    ! local
    integer(I4B) :: i, array_size
 
    call this%TspMvtType%mvt_rp()
 
    array_size = 0
    do i = 1, this%mvrbudobj%nbudterm
      array_size = array_size + this%mvrbudobj%budterm(i)%maxlist
    end do
 
    if (.not. associated(this%flow)) then
      allocate(this%flow(array_size))
      allocate(this%qty(array_size))
    end if
  
  end subroutine xmvt_rp

  subroutine xmvt_cf(this, cnew1, cnew2)
    class(TspExchangeMoverType) :: this
    real(DP), intent(in), dimension(:), contiguous, target :: cnew1
    real(DP), intent(in), dimension(:), contiguous, target :: cnew2
    ! local
    integer(I4B) :: i, n, idx
 
    call this%mvt_fill_mvrterm(cnew1, cnew2)
 
    ! add mvrterm data into synchronization array
    idx = 0
    do i = 1, size(this%mvrterm)
      do n = 1, size(this%mvrterm(i)%flow)
        this%flow(idx) = this%mvrterm(i)%flow(n)
        this%qty(idx) = this%mvrterm(i)%qty(n)
        idx = idx + 1
      end do
    end do
 
  end subroutine xmvt_cf

  subroutine xmvt_fc(this, cnew1, cnew2)
    class(TspExchangeMoverType) :: this
    real(DP), intent(in), dimension(:), contiguous, target :: cnew1
    real(DP), intent(in), dimension(:), contiguous, target :: cnew2
    ! local
    integer(I4B) :: i, n, idx
    
    ! get mvrterm data from synchronization array
    idx = 0
    do i = 1, size(this%mvrterm)
      do n = 1, size(this%mvrterm(i)%flow)
        this%mvrterm(i)%flow(n) = this%flow(idx)
        this%mvrterm(i)%qty(n) = this%qty(idx)
        idx = idx + 1
      end do
    end do
 
    call this%mvt_update_qmfrommvr()
 
  end subroutine xmvt_fc
 
  subroutine xmvt_da(this)
    class(TspExchangeMoverType) :: this
 
    if (associated(this%flow)) then
      deallocate(this%flow)
      deallocate(this%qty)
    end if
  end subroutine xmvt_da

end module TspExchangeMoverModule