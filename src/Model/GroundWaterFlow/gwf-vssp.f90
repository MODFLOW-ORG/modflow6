module VsspModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DHALF, &
                             LENFTYPE, LENPACKAGENAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_setptr, mem_checkin, get_isize
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use MatrixBaseModule

  implicit none
  private

  public :: vssp_create
  public :: VsspType

  character(len=LENFTYPE) :: ftype = 'VSSP'
  character(len=LENPACKAGENAME) :: text = '            VSSP'
  !
  type, extends(BndExtType) :: VsspType
    integer(I4B), dimension(:), pointer, contiguous :: ihc => null() !< connection type (0 vertical face, 1 horizontal)
    real(DP), dimension(:), pointer, contiguous :: cl1 => null() !< distance from cell center to seepage face
    real(DP), dimension(:), pointer, contiguous :: hwva => null() !< seepage face area
    real(DP), dimension(:), pointer, contiguous :: angldegx => null() !< face-normal angle with x axis (degrees)
    real(DP), dimension(:), pointer, contiguous :: gwfk11 => null() !< NPF horizontal hydraulic conductivity (K11), not owned
    real(DP), dimension(:), pointer, contiguous :: gwfk33 => null() !< NPF vertical hydraulic conductivity (K33), not owned
    integer(I4B), pointer :: gwfik33 => null() !< NPF flag indicating whether K33 was specified, not owned
    real(DP), dimension(:), pointer, contiguous :: gwfkrel => null() !< NPF relative permeability (KREL), not owned
    integer(I4B), dimension(:), pointer, contiguous :: vsfiunsat => null() !< VSF flag, 1 for Richards (unsaturated) nodes, not owned
  contains
    procedure :: allocate_scalars => vssp_allocate_scalars
    procedure :: allocate_arrays => vssp_allocate_arrays
    procedure :: source_options => vssp_source_options
    procedure :: bnd_cf => vssp_cf
    procedure :: bnd_fc => vssp_fc
    procedure :: bnd_fn => vssp_fn
    procedure :: bnd_ck => vssp_ck
    procedure :: bnd_da => vssp_da
    procedure :: define_listlabel
    procedure, private :: effective_k
  end type VsspType

contains

  !> @brief Create a New VSSP Package
  !<
  subroutine vssp_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(VsspType), pointer :: vsspobj

    ! allocate the object and assign values to object variables
    allocate (vsspobj)
    packobj => vsspobj

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

  end subroutine vssp_create

  subroutine vssp_source_options(this)
    class(VsspType), intent(inout) :: this

    ! source common bound options
    call this%BndExtType%source_options()

  end subroutine vssp_source_options

  !> @brief Allocate scalars
  !<
  subroutine vssp_allocate_scalars(this)
    class(VsspType) :: this !< this instance

    ! base allocate
    call this%BndExtType%allocate_scalars()

  end subroutine vssp_allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine vssp_allocate_arrays(this, nodelist, auxvar)
    class(VsspType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! local
    integer(I4B) :: isize

    ! call base type allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)

    ! set vssp input context pointers
    call mem_setptr(this%ihc, 'IHC', this%input_mempath)
    call mem_setptr(this%cl1, 'CL1', this%input_mempath)
    call mem_setptr(this%hwva, 'HWVA', this%input_mempath)
    call mem_setptr(this%angldegx, 'ANGLDEGX', this%input_mempath)

    ! checkin vssp input context pointers
    call mem_checkin(this%ihc, 'IHC', this%memoryPath, &
                     'IHC', this%input_mempath)
    call mem_checkin(this%cl1, 'CL1', this%memoryPath, &
                     'CL1', this%input_mempath)
    call mem_checkin(this%hwva, 'HWVA', this%memoryPath, &
                     'HWVA', this%input_mempath)
    call mem_checkin(this%angldegx, 'ANGLDEGX', this%memoryPath, &
                     'ANGLDEGX', this%input_mempath)

    ! set pointers to the NPF hydraulic conductivity arrays (not owned).  The
    ! LAK-style approximation reads K11/K33 directly rather than resolving the
    ! full anisotropic tensor in the seepage-face-normal direction.
    call mem_setptr(this%gwfk11, 'K11', &
                    create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfk33, 'K33', &
                    create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%gwfik33, 'IK33', &
                    create_mem_path(this%name_model, 'NPF'))

    ! set pointer to the NPF relative permeability (not owned).  KREL is 1.0
    ! for saturated cells and is reduced by the VSF unsaturated formulation.
    call mem_setptr(this%gwfkrel, 'KREL', &
                    create_mem_path(this%name_model, 'NPF'))

    ! VSSP may only be applied to Richards (VSF) cells, so the VSF package must
    ! be active.  Point at its IUNSAT flag array when present; the per-cell
    ! requirement is verified in vssp_ck.
    call get_isize('IUNSAT', create_mem_path(this%name_model, 'VSF'), isize)
    if (isize < 0) then
      write (errmsg, '(a)') &
        'The VSSP package can only be used with the VSF (Richards) package, &
        &but no VSF package is active in this model.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    else
      call mem_setptr(this%vsfiunsat, 'IUNSAT', &
                      create_mem_path(this%name_model, 'VSF'))
    end if

  end subroutine vssp_allocate_arrays

  !> @brief Effective hydraulic conductivity of cell n for the seepage face.
  !!
  !! Approximates the conductivity the way the LAK package does: the NPF grid-
  !! aligned principal conductivities are used directly, K33 (vertical) for a
  !! vertical connection (IHC = 0) and K11 (horizontal) for a horizontal
  !! connection (IHC = 1).  When K33 was not specified (IK33 = 0) the vertical
  !! value falls back to K11.  ANGLDEGX and horizontal anisotropy (K22) are not
  !! resolved by this approximation.
  !<
  function effective_k(this, n, ihc) result(hy)
    class(VsspType) :: this
    integer(I4B), intent(in) :: n !< reduced node number
    integer(I4B), intent(in) :: ihc !< connection type
    real(DP) :: hy

    if (ihc == 0) then
      ! vertical seepage face: vertical hydraulic conductivity
      if (this%gwfik33 == 0) then
        hy = this%gwfk11(n)
      else
        hy = this%gwfk33(n)
      end if
    else
      ! horizontal seepage face: horizontal hydraulic conductivity
      hy = this%gwfk11(n)
    end if
  end function effective_k

  !> @brief Formulate the HCOF and RHS terms
  !<
  subroutine vssp_cf(this)
    class(VsspType) :: this
    ! local
    integer(I4B) :: i, node
    real(DP) :: z, head, cond

    if (this%nbound .eq. 0) return

    ! Calculate hcof and rhs for each seepage face
    do i = 1, this%nbound
      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      ! seepage-face elevation datum (face centroid) used for the
      ! pressure-head switch psi = head - z
      z = DHALF * (this%dis%bot(node) + this%dis%top(node))
      head = this%xnew(node)
      if (head > z) then
        ! saturated at the face: LAK-style K11/K33 approximation scaled by the
        ! NPF relative permeability (krel) of the Richards cell
        cond = this%gwfkrel(node) * this%effective_k(node, this%ihc(i)) * &
               this%hwva(i) / this%cl1(i)
        this%hcof(i) = -cond
        this%rhs(i) = -cond * z
      else
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
      end if

    end do

  end subroutine vssp_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine vssp_fc(this, rhs, ia, idxglo, matrix_sln)
    class(VsspType) :: this
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

      ! If mover is active and this boundary is discharging,
      ! store available water (as positive value).
      ! TODO_VSF: implement mover
      if (this%imover == 1) then
        write (errmsg, '(a,a)') "Seepage Mover not supported: ", this%packName
        call store_error(errmsg, terminate=.true.)
      end if
    end do

  end subroutine vssp_fc

  !> @brief TODO_VSSP: Fill newton terms
  !<
  subroutine vssp_fn(this, rhs, ia, idxglo, matrix_sln)
    class(VsspType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    return
  end subroutine vssp_fn

  !> @brief Verify that every seepage cell is a Richards (VSF) cell
  !!
  !! The VSSP package may only be applied to unsaturated (Richards) cells, which
  !! are flagged by the VSF package through its IUNSAT array (1 for a Richards
  !! node).  A missing VSF package is reported in vssp_allocate_arrays.
  !<
  subroutine vssp_ck(this)
    class(VsspType), intent(inout) :: this
    ! local
    character(len=30) :: nodestr
    integer(I4B) :: i, node
    ! formats
    character(len=*), parameter :: fmtvssperr = &
      "('VSSP BOUNDARY ',i0,' IS APPLIED TO CELL ',a,', which is not a &
      &Richards (VSF) cell.  Set IUNSAT for this cell in the VSF package.')"

    if (.not. associated(this%vsfiunsat)) return

    do i = 1, this%nbound
      node = this%nodelist(i)
      if (node < 1) cycle
      if (this%vsfiunsat(node) /= 1) then
        call this%dis%noder_to_string(node, nodestr)
        write (errmsg, fmt=fmtvssperr) i, trim(nodestr)
        call store_error(errmsg)
      end if
    end do

    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

  end subroutine vssp_ck

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !< option is used
  subroutine define_listlabel(this)

    class(VsspType), intent(inout) :: this

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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'IHC'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'CL1'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'FACE AREA'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'ANGLDEGX'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if

  end subroutine define_listlabel

  !> @brief Deallocate memory
  !<
  subroutine vssp_da(this)
    class(VsspType) :: this

    call this%BndExtType%bnd_da()

    call mem_deallocate(this%ihc, 'IHC', this%memoryPath)
    call mem_deallocate(this%cl1, 'CL1', this%memoryPath)
    call mem_deallocate(this%hwva, 'HWVA', this%memoryPath)
    call mem_deallocate(this%angldegx, 'ANGLDEGX', this%memoryPath)

    ! NPF conductivity arrays are owned by NPF; just disassociate
    nullify (this%gwfk11)
    nullify (this%gwfk33)
    nullify (this%gwfik33)
    nullify (this%gwfkrel)

    ! VSF IUNSAT flag array is owned by VSF; just disassociate
    nullify (this%vsfiunsat)

  end subroutine vssp_da

end module VsspModule
