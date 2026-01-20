module NumericalModelModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENBUDTXT, LENPAKLOC
  use BaseModelModule, only: BaseModelType
  use BaseDisModule, only: DisBaseType
  use SparseModule, only: sparsematrix
  use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType
  use ListModule, only: ListType
  use VersionModule, only: write_listfile_header
  use MatrixBaseModule
  use VectorBaseModule

  implicit none
  private
  public :: NumericalModelType, AddNumericalModelToList, &
            GetNumericalModelFromList

  type, extends(BaseModelType) :: NumericalModelType
    character(len=LINELENGTH), pointer :: filename => null() !input file name
    integer(I4B), pointer :: neq => null() !number of equations
    integer(I4B), pointer :: nja => null() !number of connections
    integer(I4B), pointer :: moffset => null() !offset of this model in the solution
    integer(I4B), pointer :: icnvg => null() !convergence flag
    integer(I4B), dimension(:), pointer, contiguous :: ia => null() !csr row pointer
    integer(I4B), dimension(:), pointer, contiguous :: ja => null() !csr columns
    real(DP), dimension(:), pointer, contiguous :: x => null() !dependent variable (head, conc, etc)
    real(DP), dimension(:), pointer, contiguous :: rhs => null() !right-hand side vector
    real(DP), dimension(:), pointer, contiguous :: cond => null() !conductance matrix
    integer(I4B), dimension(:), pointer, contiguous :: idxglo => null() !pointer to position in solution matrix
    real(DP), dimension(:), pointer, contiguous :: xold => null() !dependent variable for previous timestep
    real(DP), dimension(:), pointer, contiguous :: flowja => null() !intercell flows
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !ibound array
    !
    ! -- Derived types
    type(ListType), pointer :: bndlist => null() !array of boundary packages for this model
    class(DisBaseType), pointer :: dis => null() !discretization object

  contains
    !
    ! -- Required for all models (override procedures defined in BaseModelType)
    procedure :: model_df
    procedure :: model_ar
    procedure :: model_fp
    procedure :: model_da
    !
    ! -- Methods specific to a numerical model
    procedure :: model_ac
    procedure :: model_mc
    procedure :: model_rp
    procedure :: model_ad
    procedure :: model_reset
    procedure :: model_solve
    procedure :: model_cf
    procedure :: model_fc
    procedure :: model_ptcchk
    procedure :: model_ptc
    procedure :: model_nr
    procedure :: model_cc
    procedure :: model_nur
    procedure :: model_cq
    procedure :: model_bd
    procedure :: model_bdcalc
    procedure :: model_bdsave
    procedure :: model_ot
    procedure :: model_bdentry
    !
    ! -- Utility methods
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: set_moffset
    procedure :: set_idsoln
    procedure :: set_xptr
    procedure :: set_rhsptr
    procedure :: set_iboundptr
    procedure :: get_mrange
    procedure :: get_mcellid
    procedure :: get_mnodeu
    procedure :: get_iasym
    procedure :: get_idv_scale
    procedure :: create_lstfile
  end type NumericalModelType

contains
  !
  ! -- Type-bound procedures for a numerical model
  !
  subroutine model_df(this)
    class(NumericalModelType) :: this
  end subroutine model_df

  subroutine model_ac(this, sparse)
    class(NumericalModelType) :: this
    type(sparsematrix), intent(inout) :: sparse
  end subroutine model_ac

  subroutine model_mc(this, matrix_sln)
    class(NumericalModelType) :: this
    class(MatrixBaseType), pointer :: matrix_sln
  end subroutine model_mc

  subroutine model_ar(this)
    class(NumericalModelType) :: this
  end subroutine model_ar

  subroutine model_rp(this)
    class(NumericalModelType) :: this
  end subroutine model_rp

  subroutine model_ad(this)
    class(NumericalModelType) :: this
  end subroutine model_ad

  subroutine model_reset(this)
    use BndModule, only: BndType, GetBndFromList
    class(NumericalModelType) :: this
    ! local
    class(BndType), pointer :: packobj
    integer(I4B) :: ip

    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_reset()
    end do

  end subroutine model_reset

  subroutine model_solve(this)
    class(NumericalModelType) :: this
  end subroutine model_solve

  subroutine model_cf(this, kiter)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: kiter
  end subroutine model_cf

  subroutine model_fc(this, kiter, matrix_sln, inwtflag)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in) :: inwtflag
  end subroutine model_fc

  subroutine model_ptcchk(this, iptc)
    class(NumericalModelType) :: this
    integer(I4B), intent(inout) :: iptc
    iptc = 0
  end subroutine model_ptcchk

  subroutine model_ptc(this, vec_residual, iptc, ptcf)
    class(NumericalModelType) :: this
    class(VectorBaseType), pointer :: vec_residual
    integer(I4B), intent(inout) :: iptc
    real(DP), intent(inout) :: ptcf
  end subroutine model_ptc

  subroutine model_nr(this, kiter, matrix, inwtflag)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix
    integer(I4B), intent(in) :: inwtflag
  end subroutine model_nr

  subroutine model_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
  end subroutine model_cc

  subroutine model_nur(this, neqmod, x, xtemp, dx, inewtonur, dxmax, locmax)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: neqmod
    real(DP), dimension(neqmod), intent(inout) :: x
    real(DP), dimension(neqmod), intent(in) :: xtemp
    real(DP), dimension(neqmod), intent(inout) :: dx
    integer(I4B), intent(inout) :: inewtonur
    real(DP), intent(inout) :: dxmax
    integer(I4B), intent(inout) :: locmax
  end subroutine model_nur

  subroutine model_cq(this, icnvg, isuppress_output)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
  end subroutine model_cq

  subroutine model_bd(this, icnvg, isuppress_output)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output
  end subroutine model_bd

  subroutine model_bdcalc(this, icnvg)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: icnvg
  end subroutine model_bdcalc

  subroutine model_bdsave(this, icnvg)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: icnvg
  end subroutine model_bdsave

  subroutine model_ot(this)
    class(NumericalModelType) :: this
  end subroutine model_ot

  subroutine model_bdentry(this, budterm, budtxt, rowlabel)
    class(NumericalModelType) :: this
    real(DP), dimension(:, :), intent(in) :: budterm
    character(len=LENBUDTXT), dimension(:), intent(in) :: budtxt
    character(len=*), intent(in) :: rowlabel
  end subroutine model_bdentry

  subroutine model_fp(this)
    class(NumericalModelType) :: this
  end subroutine model_fp

  subroutine model_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    class(NumericalModelType) :: this

    ! -- Scalars
    call mem_deallocate(this%neq)
    call mem_deallocate(this%nja)
    call mem_deallocate(this%icnvg)
    call mem_deallocate(this%moffset)
    deallocate (this%filename)
    !
    ! -- Arrays
    call mem_deallocate(this%xold)
    call mem_deallocate(this%flowja, 'FLOWJA', this%memoryPath)
    call mem_deallocate(this%idxglo)
    !
    ! -- derived types
    call this%bndlist%Clear()
    deallocate (this%bndlist)
    !
    ! -- nullify pointers
    call mem_deallocate(this%x, 'X', this%memoryPath)
    call mem_deallocate(this%rhs, 'RHS', this%memoryPath)
    call mem_deallocate(this%ibound, 'IBOUND', this%memoryPath)
    !
    ! -- BaseModelType
    call this%BaseModelType%model_da()
    !
  end subroutine model_da

  subroutine set_moffset(this, moffset)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: moffset
    this%moffset = moffset
  end subroutine set_moffset

  subroutine get_mrange(this, mstart, mend)
    class(NumericalModelType) :: this
    integer(I4B), intent(inout) :: mstart
    integer(I4B), intent(inout) :: mend
    mstart = this%moffset + 1
    mend = mstart + this%neq - 1
  end subroutine get_mrange

  subroutine set_idsoln(this, id)
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: id
    this%idsoln = id
  end subroutine set_idsoln

  subroutine allocate_scalars(this, modelname)
    use MemoryManagerModule, only: mem_allocate
    class(NumericalModelType) :: this
    character(len=*), intent(in) :: modelname
    !
    ! -- allocate basetype members
    call this%BaseModelType%allocate_scalars(modelname)
    !
    ! -- allocate members from this type
    call mem_allocate(this%neq, 'NEQ', this%memoryPath)
    call mem_allocate(this%nja, 'NJA', this%memoryPath)
    call mem_allocate(this%icnvg, 'ICNVG', this%memoryPath)
    call mem_allocate(this%moffset, 'MOFFSET', this%memoryPath)
    allocate (this%filename)
    allocate (this%bndlist)
    !
    this%filename = ''
    this%neq = 0
    this%nja = 0
    this%icnvg = 0
    this%moffset = 0
  end subroutine allocate_scalars

  subroutine allocate_arrays(this)
    use ConstantsModule, only: DZERO
    use MemoryManagerModule, only: mem_allocate
    class(NumericalModelType) :: this
    integer(I4B) :: i
    !
    call mem_allocate(this%xold, this%neq, 'XOLD', this%memoryPath)
    call mem_allocate(this%flowja, this%nja, 'FLOWJA', this%memoryPath)
    call mem_allocate(this%idxglo, this%nja, 'IDXGLO', this%memoryPath)
    !
    ! -- initialize
    do i = 1, size(this%flowja)
      this%flowja(i) = DZERO
    end do
  end subroutine allocate_arrays

  subroutine set_xptr(this, xsln, sln_offset, varNameTgt, memPathTgt)
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(NumericalModelType) :: this
    real(DP), dimension(:), pointer, contiguous, intent(in) :: xsln
    integer(I4B) :: sln_offset
    character(len=*), intent(in) :: varNameTgt
    character(len=*), intent(in) :: memPathTgt
    ! -- local
    integer(I4B) :: offset
    ! -- code
    offset = this%moffset - sln_offset
    this%x => xsln(offset + 1:offset + this%neq)
    call mem_checkin(this%x, 'X', this%memoryPath, varNameTgt, memPathTgt)
  end subroutine set_xptr

  subroutine set_rhsptr(this, rhssln, sln_offset, varNameTgt, memPathTgt)
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(NumericalModelType) :: this
    real(DP), dimension(:), pointer, contiguous, intent(in) :: rhssln
    integer(I4B) :: sln_offset
    character(len=*), intent(in) :: varNameTgt
    character(len=*), intent(in) :: memPathTgt
    ! -- local
    integer(I4B) :: offset
    ! -- code
    offset = this%moffset - sln_offset
    this%rhs => rhssln(offset + 1:offset + this%neq)
    call mem_checkin(this%rhs, 'RHS', this%memoryPath, varNameTgt, memPathTgt)
  end subroutine set_rhsptr

  subroutine set_iboundptr(this, iboundsln, sln_offset, varNameTgt, memPathTgt)
    use MemoryManagerModule, only: mem_checkin
    ! -- dummy
    class(NumericalModelType) :: this
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: iboundsln
    integer(I4B) :: sln_offset
    character(len=*), intent(in) :: varNameTgt
    character(len=*), intent(in) :: memPathTgt
    ! -- local
    integer(I4B) :: offset
    ! -- code
    offset = this%moffset - sln_offset
    this%ibound => iboundsln(offset + 1:offset + this%neq)
    call mem_checkin(this%ibound, 'IBOUND', this%memoryPath, varNameTgt, &
                     memPathTgt)
  end subroutine set_iboundptr

  subroutine get_mcellid(this, node, mcellid)
    use BndModule, only: BndType, GetBndFromList
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: node
    character(len=*), intent(inout) :: mcellid
    ! -- local
    character(len=20) :: cellid
    integer(I4B) :: ip, ipaknode, istart, istop
    class(BndType), pointer :: packobj

    if (node < 1) then
      cellid = ''
    else if (node <= this%dis%nodes) then
      call this%dis%noder_to_string(node, cellid)
    else
      cellid = '***ERROR***'
      ipaknode = node - this%dis%nodes
      istart = 1
      do ip = 1, this%bndlist%Count()
        packobj => GetBndFromList(this%bndlist, ip)
        if (packobj%npakeq == 0) cycle
        istop = istart + packobj%npakeq - 1
        if (istart <= ipaknode .and. ipaknode <= istop) then
          write (cellid, '(a, a, a, i0, a, i0, a)') '(', &
            trim(packobj%filtyp), '_', &
            packobj%ibcnum, '-', ipaknode - packobj%ioffset, ')'
          exit
        end if
        istart = istop + 1
      end do
    end if
    write (mcellid, '(i0, a, a, a, a)') this%id, '_', this%macronym, '-', &
      trim(adjustl(cellid))
  end subroutine get_mcellid

  subroutine get_mnodeu(this, node, nodeu)
    use BndModule, only: BndType, GetBndFromList
    class(NumericalModelType) :: this
    integer(I4B), intent(in) :: node
    integer(I4B), intent(inout) :: nodeu
    ! -- local
    if (node <= this%dis%nodes) then
      nodeu = this%dis%get_nodeuser(node)
    else
      nodeu = -(node - this%dis%nodes)
    end if
  end subroutine get_mnodeu

  function get_iasym(this) result(iasym)
    class(NumericalModelType) :: this
    integer(I4B) :: iasym
    iasym = 0
  end function get_iasym

  function get_idv_scale(this) result(idv_scale)
    class(NumericalModelType) :: this
    integer(I4B) :: idv_scale
    idv_scale = 0
  end function get_idv_scale

  function CastAsNumericalModelClass(obj) result(res)
    implicit none
    class(*), pointer, intent(inout) :: obj
    class(NumericalModelType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (NumericalModelType)
      res => obj
    end select
  end function CastAsNumericalModelClass

  subroutine AddNumericalModelToList(list, model)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    class(NumericalModelType), pointer, intent(inout) :: model
    ! -- local
    class(*), pointer :: obj
    !
    obj => model
    call list%Add(obj)
  end subroutine AddNumericalModelToList

  function GetNumericalModelFromList(list, idx) result(res)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    class(NumericalModelType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsNumericalModelClass(obj)
  end function GetNumericalModelFromList

  subroutine create_lstfile(this, lst_fname, model_fname, defined, headertxt)
    ! -- modules
    use KindModule, only: LGP
    use InputOutputModule, only: openfile, getunit
    ! -- dummy
    class(NumericalModelType) :: this
    character(len=*), intent(inout) :: lst_fname
    character(len=*), intent(in) :: model_fname
    logical(LGP), intent(in) :: defined
    character(len=*), intent(in) :: headertxt
    ! -- local
    integer(I4B) :: i, istart, istop
    !
    ! -- set list file name if not provided
    if (.not. defined) then
      !
      ! -- initialize
      lst_fname = ' '
      istart = 0
      istop = len_trim(model_fname)
      !
      ! -- identify '.' character position from back of string
      do i = istop, 1, -1
        if (model_fname(i:i) == '.') then
          istart = i
          exit
        end if
      end do
      !
      ! -- if not found start from string end
      if (istart == 0) istart = istop + 1
      !
      ! -- set list file name
      lst_fname = model_fname(1:istart)
      istop = istart + 3
      lst_fname(istart:istop) = '.lst'
    end if
    !
    ! -- create the list file
    this%iout = getunit()
    call openfile(this%iout, 0, lst_fname, 'LIST', filstat_opt='REPLACE')
    !
    ! -- write list file header
    call write_listfile_header(this%iout, headertxt)
  end subroutine create_lstfile

end module NumericalModelModule
