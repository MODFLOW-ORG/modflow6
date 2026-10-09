module SgdModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LENFTYPE, LENPACKAGENAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use MemoryManagerModule, only: mem_setptr, mem_checkin, mem_deallocate
  use MemoryHelperModule, only: create_mem_path
  use ConductanceProviderModule, only: ConductanceProviderType
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use MatrixBaseModule

  implicit none
  private

  public :: sgd_create
  public :: SgdType

  character(len=LENFTYPE) :: ftype = 'SGD'
  character(len=LENPACKAGENAME) :: text = '             SGD'
  !
  !> @brief Specified Gradient (SGD) boundary package.
  !!
  !! Imposes a user-specified hydraulic gradient across boundary faces,
  !! using a directional effective-K provider to compute the resulting flows.
  !<
  type, extends(BndExtType) :: SgdType
    real(DP), dimension(:), pointer, contiguous :: gradx => null() !< specified gradient x component (flow direction)
    real(DP), dimension(:), pointer, contiguous :: grady => null() !< specified gradient y component (flow direction)
    real(DP), dimension(:), pointer, contiguous :: gradz => null() !< specified gradient z component (flow direction)
    real(DP), dimension(:), pointer, contiguous :: hwva => null() !< boundary face area
    class(ConductanceProviderType), pointer :: cond_provider => null() !< directional effective-K provider (injected by the model)
  contains
    procedure :: allocate_arrays => sgd_allocate_arrays
    procedure :: bnd_cf => sgd_cf
    procedure :: bnd_fc => sgd_fc
    procedure :: bnd_fn => sgd_fn
    procedure :: bnd_da => sgd_da
    procedure :: define_listlabel
    procedure, private :: cond_factor
  end type SgdType

contains

  !> @brief Create a New SGD Package
  !<
  subroutine sgd_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(SgdType), pointer :: sgdobj

    ! allocate the object and assign values to object variables
    allocate (sgdobj)
    packobj => sgdobj

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

  end subroutine sgd_create

  !> @brief Allocate arrays
  !<
  subroutine sgd_allocate_arrays(this, nodelist, auxvar)
    class(SgdType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar

    ! call base type allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)

    ! set sgd input context pointers
    call mem_setptr(this%gradx, 'GRADX', this%input_mempath)
    call mem_setptr(this%grady, 'GRADY', this%input_mempath)
    call mem_setptr(this%gradz, 'GRADZ', this%input_mempath)
    call mem_setptr(this%hwva, 'HWVA', this%input_mempath)

    ! checkin sgd input context pointers
    call mem_checkin(this%gradx, 'GRADX', this%memoryPath, &
                     'GRADX', this%input_mempath)
    call mem_checkin(this%grady, 'GRADY', this%memoryPath, &
                     'GRADY', this%input_mempath)
    call mem_checkin(this%gradz, 'GRADZ', this%memoryPath, &
                     'GRADZ', this%input_mempath)
    call mem_checkin(this%hwva, 'HWVA', this%memoryPath, &
                     'HWVA', this%input_mempath)

  end subroutine sgd_allocate_arrays

  !> @brief Head-independent conductance factor Keff*|grad|*A of boundary i.
  !!
  !! Keff is the saturated directional effective hydraulic conductivity along
  !! the unit gradient, resolved through the model conductance provider. The
  !! relative permeability (saturation) factor is applied separately so this
  !< factor can be reused by both the formulate and Newton routines.
  function cond_factor(this, i, node) result(cond)
    class(SgdType) :: this
    integer(I4B), intent(in) :: i !< boundary number
    integer(I4B), intent(in) :: node !< reduced node number
    real(DP) :: cond
    ! local
    integer(I4B) :: ihc
    real(DP) :: gmag, keff
    real(DP), dimension(3) :: ghat

    cond = DZERO
    gmag = sqrt(this%gradx(i)**2 + this%grady(i)**2 + this%gradz(i)**2)
    if (gmag <= DZERO) return

    ! unit gradient (flow direction) and vertical/horizontal anisotropy branch
    ghat = [this%gradx(i), this%grady(i), this%gradz(i)] / gmag
    if (this%gradx(i) == DZERO .and. this%grady(i) == DZERO) then
      ihc = 0
    else
      ihc = 1
    end if

    keff = this%cond_provider%eff_hy(node, ihc, ghat)
    cond = keff * gmag * this%hwva(i)
  end function cond_factor

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! A specified gradient drives a Darcy flux q = krel*Keff*|grad| that leaves
  !! the model across a boundary face of area HWVA. The flux is head
  !! independent (a Neumann term) apart from the lagged relative permeability,
  !< so it contributes to RHS only. Free drainage is grad = (0, 0, -1).
  subroutine sgd_cf(this)
    class(SgdType) :: this
    ! local
    integer(I4B) :: i, node
    real(DP) :: kr, cond

    if (this%nbound == 0) return

    do i = 1, this%nbound
      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      cond = this%cond_factor(i, node)
      kr = this%cond_provider%krel(node)

      this%hcof(i) = DZERO
      this%rhs(i) = kr * cond
    end do

  end subroutine sgd_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine sgd_fc(this, rhs, ia, idxglo, matrix_sln)
    class(SgdType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, n, ipos

    ! no mover support
    if (this%imover == 1) then
      write (errmsg, '(a,a)') "SGD Mover not supported: ", this%packName
      call store_error(errmsg, terminate=.true.)
    end if

    ! Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
    end do

  end subroutine sgd_fc

  !> @brief Add Newton-Raphson terms for the package into the solution.
  !!
  !! The boundary flux depends on head only through the relative permeability
  !! krel(h), so the RHS term is q(h) = krel(h) * cond with cond head
  !! independent. The Newton contribution is the linearization of that term:
  !! the diagonal gains -d(q)/dh and the RHS gains -d(q)/dh * h, where
  !< d(q)/dh = cond * dkrel/dh is supplied by the conductance provider.
  subroutine sgd_fn(this, rhs, ia, idxglo, matrix_sln)
    class(SgdType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, node, ipos
    real(DP) :: dkrdh, cond, drterm

    do i = 1, this%nbound
      node = this%nodelist(i)
      if (this%ibound(node) <= 0) cycle

      ! derivative of relative permeability (0 for saturated cells -> no term)
      dkrdh = this%cond_provider%dkrel_dh(node)
      if (dkrdh == DZERO) cycle

      cond = this%cond_factor(i, node)
      drterm = cond * dkrdh

      ipos = ia(node)
      call matrix_sln%add_value_pos(idxglo(ipos), -drterm)
      rhs(node) = rhs(node) - drterm * this%xnew(node)
    end do

  end subroutine sgd_fn

  !> @brief Deallocate memory
  !<
  subroutine sgd_da(this)
    class(SgdType) :: this

    ! deallocate base
    call this%BndExtType%bnd_da()

    ! deallocate input context pointers
    call mem_deallocate(this%gradx, 'GRADX', this%memoryPath)
    call mem_deallocate(this%grady, 'GRADY', this%memoryPath)
    call mem_deallocate(this%gradz, 'GRADZ', this%memoryPath)
    call mem_deallocate(this%hwva, 'HWVA', this%memoryPath)

    ! provider is owned by the model; just disassociate
    nullify (this%cond_provider)

  end subroutine sgd_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !< option is used
  subroutine define_listlabel(this)
    class(SgdType), intent(inout) :: this

    ! create the header list label
    this%listlabel = trim(this%filtyp)//' NO.'
    if (this%dis%ndim == 3) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'ROW'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'COL'
    elseif (this%dis%ndim == 2) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'CELL2D'
    else
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'NODE'
    end if
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'GRADX'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'GRADY'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'GRADZ'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'HWVA'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if

  end subroutine define_listlabel

end module SgdModule
