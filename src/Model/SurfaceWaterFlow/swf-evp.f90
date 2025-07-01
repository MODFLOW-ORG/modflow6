!> @brief This module contains the evaporation (EVP) package methods
!!
!! This module can be used to represent evaporation onto streams and
!! overland flow cells.
!<
module SwfEvpModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, LENFTYPE, LENPACKAGENAME, MAXCHARLEN, &
                             LINELENGTH, DONE, DHALF, DEM6, DPREC
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use SimModule, only: store_error, store_error_filename, count_errors
  use SimVariablesModule, only: errmsg
  use ObsModule, only: DefaultObsIdProcessor
  use TimeArraySeriesLinkModule, only: TimeArraySeriesLinkType
  use BlockParserModule, only: BlockParserType
  use CharacterStringModule, only: CharacterStringType
  use SmoothingModule, only: sQSaturation
  use MathUtilModule, only: get_perturbation
  use MatrixBaseModule
  use GeomUtilModule, only: get_node
  use BaseDisModule, only: DisBaseType
  use Disv1dModule, only: Disv1dType
  use SwfDfwModule, only: SwfDfwType
  use SwfCxsModule, only: SwfCxsType

  implicit none

  private
  public :: evp_create

  character(len=LENFTYPE) :: ftype = 'EVP'
  character(len=LENPACKAGENAME) :: text = '             EVP'
  ! character(len=LENPACKAGENAME) :: texta = '            EVPA'

  type, extends(BndExtType) :: SwfEvpType
    real(DP), dimension(:), pointer, contiguous :: evaporation => null() !< boundary evaporation array
    integer(I4B), pointer :: iflowred => null() !< flag that indicates evaporation will be shut off when depth is less than reduction depth
    real(DP), pointer :: reduction_depth => null() !< depth below which evaporation is reduced
    logical, pointer, private :: read_as_arrays

    ! pointers to other objects
    type(SwfDfwType), pointer :: dfw
    type(SwfCxsType), pointer :: cxs

  contains

    procedure :: evp_allocate_scalars
    procedure :: allocate_arrays => evp_allocate_arrays
    procedure :: source_options => evp_source_options
    procedure :: source_dimensions => evp_source_dimensions
    procedure :: log_evp_options
    procedure :: read_initial_attr => evp_read_initial_attr
    procedure :: bnd_rp => evp_rp
    procedure :: bnd_ck => evp_ck
    procedure :: bnd_cf => evp_cf
    procedure :: bnd_fc => evp_fc
    procedure :: bnd_da => evp_da
    procedure :: define_listlabel => evp_define_listlabel
    procedure :: bound_value => evp_bound_value
    procedure :: default_nodelist
    procedure, private :: reach_length_pointer
    procedure, private :: get_qevp
    procedure, private :: get_evap_reduce_mult
    ! for observations
    procedure, public :: bnd_obs_supported => evp_obs_supported
    procedure, public :: bnd_df_obs => evp_df_obs

  end type SwfEvpType

contains

  !> @brief Create a Evaporation Package
  !<
  subroutine evp_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath, dis, dfw, cxs)
    ! dummy
    class(BndType), pointer :: packobj !< pointer to default package type
    integer(I4B), intent(in) :: id !< package id
    integer(I4B), intent(in) :: ibcnum !< boundary condition number
    integer(I4B), intent(in) :: inunit !< unit number of CDB package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< model name
    character(len=*), intent(in) :: pakname !< package name
    character(len=*), intent(in) :: mempath !< input mempath
    class(DisBaseType), pointer, intent(inout) :: dis !< the pointer to the discretization
    type(SwfDfwType), pointer, intent(in) :: dfw !< the pointer to the dfw package
    type(SwfCxsType), pointer, intent(in) :: cxs !< the pointer to the cxs package
    ! local
    type(SwfEvpType), pointer :: evpobj

    ! allocate evaporation object and scalar variables
    allocate (evpobj)
    packobj => evpobj

    ! create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text

    ! allocate scalars
    call evpobj%evp_allocate_scalars()

    ! initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'DFW')

    ! store pointer to dis
    evpobj%dis => dis

    ! store pointer to dfw
    evpobj%dfw => dfw

    ! store pointer to cxs
    evpobj%cxs => cxs
  end subroutine evp_create

  !> @brief Allocate scalar members
  !<
  subroutine evp_allocate_scalars(this)
    ! dummy
    class(SwfEvpType), intent(inout) :: this

    ! allocate base scalars
    call this%BndExtType%allocate_scalars()

    ! allocate internal members
    call mem_allocate(this%iflowred, 'IFLOWRED', this%memoryPath)
    call mem_allocate(this%reduction_depth, 'REDUCTION_DEPTH', this%memoryPath)
    allocate (this%read_as_arrays)

    ! Set values
    this%iflowred = 1
    this%reduction_depth = DEM6
    this%read_as_arrays = .false.
  end subroutine evp_allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine evp_allocate_arrays(this, nodelist, auxvar)
    ! modules
    use MemoryManagerModule, only: mem_setptr, mem_checkin
    ! dummy
    class(SwfEvpType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar

    ! allocate base arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)

    ! set input context pointers
    call mem_setptr(this%evaporation, 'EVAPORATION', this%input_mempath)

    ! checkin input context pointers
    call mem_checkin(this%evaporation, 'EVAPORATION', this%memoryPath, &
                     'EVAPORATION', this%input_mempath)
  end subroutine evp_allocate_arrays

  !> @brief Source options specific to EVPType
  !<
  subroutine evp_source_options(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    implicit none
    ! dummy
    class(SwfEvpType), intent(inout) :: this
    ! local
    logical(LGP) :: found_readasarrays = .false.

    ! source common bound options
    call this%BndExtType%source_options()

    ! update defaults with idm sourced values
    call mem_set_value(this%read_as_arrays, 'READASARRAYS', this%input_mempath, &
                       found_readasarrays)

    ! log evp params
    call this%log_evp_options(found_readasarrays)
  end subroutine evp_source_options

  !> @brief Log options specific to SwfEvpType
  !<
  subroutine log_evp_options(this, found_readasarrays)
    implicit none
    ! dummy
    class(SwfEvpType), intent(inout) :: this
    logical(LGP), intent(in) :: found_readasarrays
    ! formats
    character(len=*), parameter :: fmtreadasarrays = &
      &"(4x, 'EVAPORATION INPUT WILL BE READ AS ARRAY(S).')"

    ! log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'

    if (found_readasarrays) then
      write (this%iout, fmtreadasarrays)
    end if

    ! close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine log_evp_options

  !> @brief Source the dimensions for this package
  !<
  subroutine evp_source_dimensions(this)
    ! dummy
    class(SwfEvpType), intent(inout) :: this

    if (this%read_as_arrays) then

      ! Set maxbound to the number of cells per layer, which is simply
      ! nrow * ncol for a dis2d grid, and nodesuser for disv2d and disv1d
      this%maxbound = this%dis%get_ncpl()

      ! verify dimensions were set
      if (this%maxbound <= 0) then
        write (errmsg, '(a)') &
          'MAXBOUND must be an integer greater than zero.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if

    else

      ! source maxbound
      call this%BndExtType%source_dimensions()

    end if

    ! Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
  end subroutine evp_source_dimensions

  !> @brief Part of allocate and read
  !<
  subroutine evp_read_initial_attr(this)
    ! dummy
    class(SwfEvpType), intent(inout) :: this

    if (this%read_as_arrays) then
      call this%default_nodelist()
    end if
  end subroutine evp_read_initial_attr

  !> @brief Read and Prepare
  !!
  !! Read itmp and read new boundaries if itmp > 0
  !<
  subroutine evp_rp(this)
    ! modules
    use TdisModule, only: kper
    implicit none
    ! dummy
    class(SwfEvpType), intent(inout) :: this

    if (this%iper /= kper) return

    if (this%read_as_arrays) then
      ! no need to do anything because this%evaporation points directly to
      ! the input context evaporation, which is automatically updated by idm
    else
      call this%BndExtType%bnd_rp()
    end if

    ! Write the list to iout if requested
    if (this%iprpak /= 0) then
      call this%write_list()
    end if
  end subroutine evp_rp

  !> @brief Ensure evaporation is positive
  !<
  subroutine evp_ck(this)
    ! dummy
    class(SwfEvpType), intent(inout) :: this
    ! local
    character(len=30) :: nodestr
    integer(I4B) :: i, nr
    character(len=*), parameter :: fmterr = &
      &"('Specified stress ',i0, &
      &' evaporation (',g0,') is less than zero for cell', a)"

    ! Ensure evaporation rates are positive
    do i = 1, this%nbound
      nr = this%nodelist(i)
      if (nr <= 0) cycle
      if (this%evaporation(i) < DZERO) then
        call this%dis%noder_to_string(nr, nodestr)
        write (errmsg, fmt=fmterr) i, this%evaporation(i), trim(nodestr)
        call store_error(errmsg)
      end if
    end do

    ! write summary of package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine evp_ck

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Skip if no evaporation. Otherwise, calculate hcof and rhs
  !<
  subroutine evp_cf(this)
    ! dummy
    class(SwfEvpType) :: this
    ! local
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: inwt
    real(DP) :: q
    real(DP) :: qeps
    real(DP) :: eps
    real(DP) :: derv
    real(DP) :: evap
    real(DP) :: rlen
    real(DP), dimension(:), pointer :: reach_length

    ! Return if no evaporation
    if (this%nbound == 0) return

    ! Set pointer to reach_length for 1d
    reach_length => this%reach_length_pointer()
    rlen = DZERO

    ! Calculate hcof and rhs for each evaporation entry
    do i = 1, this%nbound

      ! Find the node number
      node = this%nodelist(i)

      ! cycle if nonexistent bound
      if (node <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      ! cycle if dry or constant head
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      ! Initialize hcof
      this%hcof(i) = DZERO

      ! assign evap rate in length per time and multiply by auxvar
      evap = this%evaporation(i)
      if (this%iauxmultcol > 0) then
        evap = evap * this%auxvar(this%iauxmultcol, i)
      end if

      ! get reach length for 1d channel
      if (this%dis%is_1d()) then
        rlen = reach_length(node)
      end if

      ! Calculate volumetric evaporation flow in L^3/Td and add to rhs
      q = -this%get_qevp(node, rlen, this%xnew(node), this%xold(node), evap)
      this%rhs(i) = -q

      ! Code for adding newton terms
      inwt = 1
      if (inwt == 1) then

        ! calculate perturbed q
        eps = get_perturbation(this%xnew(node))
        qeps = -this%get_qevp(node, rlen, this%xnew(node) + eps, &
                              this%xold(node), evap)

        ! calculate derivative
        derv = (qeps - q) / eps

        ! add derivative to hcof and update rhs with derivate contribution
        this%hcof(i) = derv
        this%rhs(i) = this%rhs(i) + derv * this%xnew(node)
      end if

    end do
  end subroutine evp_cf

  !> @brief Calculate qevp
  !!
  !! Calculate qevp for both channel and overland flow grids.
  !! Approximate the average water surface width of the channel
  !! as wavg = delta A over delta h, and then multiply wavg
  !! by reach length to come up with surface water area for the
  !! channel.  Reduce evaporation when depths are small and shut
  !! it off when there is no water in the cell.
  !<
  function get_qevp(this, node, rlen, snew, sold, evaporation) result(qevp)
    ! dummy
    class(SwfEvpType) :: this !< this instance
    integer(I4B), intent(in) :: node !< reduced node number
    real(DP), intent(in) :: rlen !< length of reach
    real(DP), intent(in) :: snew !< current stage in reach
    real(DP), intent(in) :: sold !< previous stage in reach
    real(DP), intent(in) :: evaporation !< evaporation rate in length per time
    ! return
    real(DP) :: qevp
    ! local
    integer(I4B) :: idcxs
    real(DP) :: depth
    real(DP) :: bt
    real(DP) :: area
    real(DP) :: anew
    real(DP) :: aold
    real(DP) :: denom
    real(DP) :: wavg
    real(DP) :: width_channel
    real(DP) :: dummy
    real(DP) :: qmult

    ! calculate depth of water
    bt = this%dis%bot(node)

    ! Determine the water surface area
    if (this%dis%is_2d()) then
      ! overland flow case
      area = this%dis%get_area(node)
    else if (this%dis%is_1d()) then
      ! channel case
      idcxs = this%dfw%idcxs(node)
      call this%dis%get_flow_width(node, node, 0, width_channel, dummy)

      depth = snew - bt
      anew = this%cxs%get_area(idcxs, width_channel, depth)
      depth = sold - bt
      aold = this%cxs%get_area(idcxs, width_channel, depth)
      wavg = this%cxs%get_wetted_top_width(idcxs, width_channel, depth)
      denom = snew - sold
      if (abs(denom) > DPREC) then
        wavg = (anew - aold) / (snew - sold)
      end if
      area = rlen * wavg

    end if

    ! Reduce the evap rate as cell goes dry
    qmult = this%get_evap_reduce_mult(snew, bt)

    ! calculate volumetric evaporation flow in L^3/T
    qevp = evaporation * area * qmult

  end function get_qevp

  !> @brief Calculate multiplier to reduce evap as depth goes to zero
  !<
  function get_evap_reduce_mult(this, stage, bottom) result(qmult)
    ! dummy
    class(SwfEvpType) :: this
    real(DP), intent(in) :: stage
    real(DP), intent(in) :: bottom
    ! return
    real(DP) :: qmult
    ! local
    real(DP) :: tp

    qmult = DONE
    if (this%iflowred == 1) then
      tp = bottom + this%reduction_depth
      qmult = sQSaturation(tp, bottom, stage)
    end if

  end function get_evap_reduce_mult

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine evp_fc(this, rhs, ia, idxglo, matrix_sln)
    ! dummy
    class(SwfEvpType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, n, ipos

    ! Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      if (n <= 0) cycle
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
    end do
  end subroutine evp_fc

  !> @brief Deallocate memory
  !<
  subroutine evp_da(this)
    ! modules
    ! dummy
    class(SwfEvpType) :: this

    ! Deallocate parent package
    call this%BndExtType%bnd_da()

    ! scalars
    call mem_deallocate(this%iflowred)
    call mem_deallocate(this%reduction_depth)
    deallocate (this%read_as_arrays)

    ! arrays
    call mem_deallocate(this%evaporation, 'EVAPORATION', this%memoryPath)

    ! pointers
    nullify (this%dis)
    nullify (this%dfw)
    nullify (this%cxs)
  end subroutine evp_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used.
  !<
  subroutine evp_define_listlabel(this)
    ! dummy
    class(SwfEvpType), intent(inout) :: this
    !
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'EVAPORATION'
!    if(this%multindex > 0) &
!      write(this%listlabel, '(a, a16)') trim(this%listlabel), 'MULTIPLIER'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine evp_define_listlabel

  !> @brief Assign default nodelist when READASARRAYS is specified.
  !<
  subroutine default_nodelist(this)
    ! dummy
    class(SwfEvpType) :: this
    ! local
    integer(I4B) :: nodeu, noder

    ! This is only called for readasarrays, so nodelist will be the size of
    ! the user grid, and will have a value of 0 for any entries where idomain
    ! is not 1
    do nodeu = 1, this%maxbound
      noder = this%dis%get_nodenumber(nodeu, 0)
      this%nodelist(nodeu) = noder
    end do

    ! Assign nbound
    this%nbound = this%maxbound

  end subroutine default_nodelist

  ! Procedures related to observations

  !> @brief
  !!
  !! Overrides BndType%bnd_obs_supported()
  !<
  logical function evp_obs_supported(this)
    implicit none
    ! dummy
    class(SwfEvpType) :: this
    evp_obs_supported = .true.
  end function evp_obs_supported

  !> @brief Implements bnd_df_obs
  !!
  !! Store observation type supported by EVP package. Overrides
  !! BndType%bnd_df_obs
  !<
  subroutine evp_df_obs(this)
    implicit none
    ! dummy
    class(SwfEvpType) :: this
    ! local
    integer(I4B) :: indx

    call this%obs%StoreObsType('evp', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine evp_df_obs

  !> @brief Return requested boundary value
  !<
  function evp_bound_value(this, col, row) result(bndval)
    ! modules
    use ConstantsModule, only: DZERO
    ! dummy
    class(SwfEvpType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! result
    real(DP) :: bndval

    select case (col)
    case (1)
      if (this%iauxmultcol > 0) then
        bndval = this%evaporation(row) * this%auxvar(this%iauxmultcol, row)
      else
        bndval = this%evaporation(row)
      end if
    case default
      errmsg = 'Programming error. EVP bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function evp_bound_value

  function reach_length_pointer(this) result(ptr)
    ! dummy
    class(SwfEvpType) :: this !< this instance
    ! return
    real(DP), dimension(:), pointer :: ptr
    ! local
    class(DisBaseType), pointer :: dis

    ptr => null()
    dis => this%dis
    select type (dis)
    type is (Disv1dType)
      ptr => dis%length
    end select

  end function reach_length_pointer

end module SwfEvpModule

