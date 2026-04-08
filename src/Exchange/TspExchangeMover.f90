module TspExchangeMoverModule
  use TspMvtModule

  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: LENVARNAME, DNODATA
  use VirtualModelModule, only: VirtualModelType, get_virtual_model
  use TspFmiModule, only: TspFmiType
  use GwtModule, only: GwtModelType
  implicit none
  private

  public :: xmvt_cr

  type, public, extends(TspMvtType) :: TspExchangeMoverType
    class(VirtualModelType), pointer :: model1 => null() !< virtual model 1
    class(VirtualModelType), pointer :: model2 => null() !< virtual model 2
    real(DP), dimension(:), pointer, contiguous :: qty_m1 => null() !< concentration, temp, ...
    real(DP), dimension(:), pointer, contiguous :: qty_m2 => null() !< concentration, temp, ...
  contains
    procedure :: mvt_rp => xmvt_rp
    procedure :: xmvt_cf
    procedure :: mvt_fc => xmvt_fc
    procedure :: mvt_da => xmvt_da
  end type TspExchangeMoverType

contains
  
  subroutine xmvt_cr(mvt, name_exg, gwt_model1, gwt_model2, &
                     gwfmodelname1, gwfmodelname2, inunit, iout)
    type(TspExchangeMoverType), pointer :: mvt    
    character(len=*), intent(in) :: name_exg !< name of the exchange
    class(GwtModelType), pointer :: gwt_model1 !< gwt model 1, can be null()
    class(GwtModelType), pointer :: gwt_model2 !< gwt model 2, can be null()
    character(len=*), intent(in) :: gwfmodelname1
    character(len=*), intent(in) :: gwfmodelname2
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    ! local
    type(TspFmiType), pointer :: fmi1, fmi2
    real(DP), pointer :: eqnsclfac !< governing equation scale factor    
    character(len=LENVARNAME) :: depvartype !< dependent variable type ('concentration' or 'temperature')


    allocate (mvt)

    fmi1 => null()
    if (associated(gwt_model1)) then
      fmi1 => gwt_model1%fmi
      eqnsclfac => gwt_model1%eqnsclfac
      depvartype = gwt_model1%depvartype
    end if
    fmi2 => null()
    if (associated(gwt_model2)) then
      fmi2 => gwt_model2%fmi
      eqnsclfac => gwt_model2%eqnsclfac
      depvartype = gwt_model2%depvartype
    end if

    mvt%model1 => get_virtual_model(gwfmodelname1)
    mvt%model2 => get_virtual_model(gwfmodelname2)

    call mvt%mvt_init(name_exg, inunit, iout, fmi1, &
                      eqnsclfac, depvartype, gwfmodelname1, &
                      gwfmodelname2, fmi2)

  end subroutine xmvt_cr

  !> @brief Read and prepare mover transport object  
  !<
  subroutine xmvt_rp(this)
    class(TspExchangeMoverType) :: this
    ! local
    integer(I4B) :: i, max_array_size
 
    call this%TspMvtType%mvt_rp()
 
    max_array_size = 0
    do i = 1, this%mvrbudobj%nbudterm
      max_array_size = max_array_size + this%mvrbudobj%budterm(i)%maxlist
    end do
 
    if (associated(this%qty_m1)) deallocate(this%qty_m1)
    allocate(this%qty_m1(max_array_size))
    if (associated(this%qty_m2)) deallocate(this%qty_m2)
    allocate(this%qty_m2(max_array_size))

    this%qty_m1 = DNODATA
    this%qty_m2 = DNODATA
  
  end subroutine xmvt_rp

  subroutine xmvt_cf(this, cnew1, cnew2)
    class(TspExchangeMoverType) :: this
    real(DP), intent(in), dimension(:), contiguous, target :: cnew1
    real(DP), intent(in), dimension(:), contiguous, target :: cnew2
    ! local
    integer(I4B) :: i, n, idx
    real(DP), dimension(:), pointer, contiguous :: prov_array
    logical(LGP) :: prov_is_local

    call this%mvt_fill_mvrterm(cnew1, cnew2)
 
    ! add mvrterm data into synchronization array
    idx = 1
    do i = 1, size(this%mvrterm)
      
      ! point to provider array
      if (this%mvrterm(i)%prov_is_m1) then
        prov_is_local = this%model1%is_local
        prov_array => this%qty_m1
      else
        prov_is_local = this%model2%is_local
        prov_array => this%qty_m2
      end if

      if (prov_is_local) then
        do n = 1, size(this%mvrterm(i)%qty)
          prov_array(idx) = this%mvrterm(i)%qty(n)
          idx = idx + 1
        end do
      end if
    end do
 
  end subroutine xmvt_cf

  subroutine xmvt_fc(this, cnew1, cnew2)
    use SimVariablesModule, only: proc_id
    class(TspExchangeMoverType) :: this
    real(DP), intent(in), dimension(:), contiguous, target :: cnew1
    real(DP), intent(in), dimension(:), contiguous, target :: cnew2
    ! local
    integer(I4B) :: i, n, idx
    real(DP), dimension(:), pointer, contiguous :: prov_array ! after synchronization when parallel
    logical(LGP) :: prov_is_remote

    ! get mvrterm data from synchronization array
    idx = 1
    do i = 1, size(this%mvrterm)

      ! point to provider array
      if (this%mvrterm(i)%prov_is_m1) then
        prov_is_remote = .not. this%model1%is_local
        prov_array => this%qty_m1
      else
        prov_is_remote = .not. this%model2%is_local
        prov_array => this%qty_m2
      end if

      if (prov_is_remote) then
        ! remote provider, so we copy it in here after sync'ing
        do n = 1, size(this%mvrterm(i)%qty)
          write(*,*) proc_id, ": moving ", prov_array(idx), this%mvrbudobj%budterm(i)%flow(n)
          this%mvrterm(i)%qty(n) = prov_array(idx)
          idx = idx + 1
        end do
      end if
    end do
 
    call this%mvt_update_qmfrommvr()
 
  end subroutine xmvt_fc
 
  subroutine xmvt_da(this)
    class(TspExchangeMoverType) :: this

    call this%TspMvtType%mvt_da()
 
    if (associated(this%qty_m1)) then
      deallocate(this%qty_m1)
    end if
  end subroutine xmvt_da

end module TspExchangeMoverModule