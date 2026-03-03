module SfbModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DHALF, DEM6, LENFTYPE, LENPACKAGENAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_setptr, mem_checkin
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use MatrixBaseModule

  implicit none
  private

  public :: sfb_create
  public :: SfbType

  character(len=LENFTYPE) :: ftype = 'SFB'
  character(len=LENPACKAGENAME) :: text = '             SFB'
  !
  type, extends(BndExtType) :: SfbType
  contains
    procedure :: bnd_cf => sfb_cf
    procedure :: bnd_fc => sfb_fc
  end type SfbType

contains

  !> @brief Create a New SFB Package
  !<
  subroutine sfb_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath)
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    ! local
    type(SfbType), pointer :: sfbobj

    ! allocate the object and assign values to object variables
    allocate (sfbobj)
    packobj => sfbobj

    ! create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text

    ! allocate scalars
    call packobj%allocate_scalars()

    ! initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  
  end subroutine sfb_create

  !> @brief Formulate the HCOF and RHS terms
  !<
  subroutine sfb_cf(this)
    class(SfbType) :: this
    ! local
    integer(I4B) :: i, node
    real(DP), dimension(:), pointer, contiguous :: npf_k33 => null()
    real(DP), dimension(:), pointer, contiguous :: npf_krel => null()

    call mem_setptr(npf_k33, 'K33', this%ictMemPath)
    call mem_setptr(npf_krel, 'KREL', this%ictMemPath)

    ! Calculate hcof and rhs for each seepage face
    do i = 1, this%nbound
      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      this%hcof(i) = DZERO
      this%rhs(i) = npf_krel(node) * npf_k33(node) * this%dis%area(node) ! times the unit gradient
      
    end do

  end subroutine sfb_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine sfb_fc(this, rhs, ia, idxglo, matrix_sln)
    class(SfbType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, n, ipos

    ! pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if

    ! Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))

      ! no mover support here
      if (this%imover == 1) then
        write (errmsg, '(a,a)') "SFB Mover not supported: ", this%packName
        call store_error(errmsg, terminate=.true.)
      end if
    end do

  end subroutine sfb_fc

end module SfbModule
