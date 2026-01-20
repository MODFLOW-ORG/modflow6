!> @brief This module contains the precipitation (PCP) package methods
!!
!! This module can be used to represent precipitation onto streams and
!! overland flow cells.
!<
module SwfPcpModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, LENFTYPE, LENPACKAGENAME, MAXCHARLEN, &
                             LINELENGTH
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use SimModule, only: store_error, store_error_filename, count_errors
  use SimVariablesModule, only: errmsg
  use ObsModule, only: DefaultObsIdProcessor
  use TimeArraySeriesLinkModule, only: TimeArraySeriesLinkType
  use BlockParserModule, only: BlockParserType
  use CharacterStringModule, only: CharacterStringType
  use MatrixBaseModule
  use GeomUtilModule, only: get_node
  use BaseDisModule, only: DisBaseType
  use Disv1dModule, only: Disv1dType
  use SwfDfwModule, only: SwfDfwType
  use SwfCxsModule, only: SwfCxsType

  implicit none

  private
  public :: pcp_create

  character(len=LENFTYPE) :: ftype = 'PCP'
  character(len=LENPACKAGENAME) :: text = '             PCP'
  ! character(len=LENPACKAGENAME) :: texta = '            PCPA'

  type, extends(BndExtType) :: SwfPcpType
    real(DP), dimension(:), pointer, contiguous :: precipitation => null() !< boundary precipitation array
    logical, pointer, private :: read_as_arrays

    ! pointers to other objects
    type(SwfDfwType), pointer :: dfw
    type(SwfCxsType), pointer :: cxs

  contains

    procedure :: pcp_allocate_scalars
    procedure :: allocate_arrays => pcp_allocate_arrays
    procedure :: source_options => pcp_source_options
    procedure :: source_dimensions => pcp_source_dimensions
    procedure :: log_pcp_options
    procedure :: read_initial_attr => pcp_read_initial_attr
    procedure :: bnd_rp => pcp_rp
    procedure :: bnd_ck => pcp_ck
    procedure :: bnd_cf => pcp_cf
    procedure :: bnd_fc => pcp_fc
    procedure :: bnd_da => pcp_da
    procedure :: define_listlabel => pcp_define_listlabel
    procedure :: bound_value => pcp_bound_value
    procedure :: default_nodelist
    procedure, private :: reach_length_pointer
    ! for observations
    procedure, public :: bnd_obs_supported => pcp_obs_supported
    procedure, public :: bnd_df_obs => pcp_df_obs

  end type SwfPcpType

contains

  !> @brief Create a Precipitation Package
  !<
  subroutine pcp_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(SwfPcpType), pointer :: pcpobj

    ! allocate precipitation object and scalar variables
    allocate (pcpobj)
    packobj => pcpobj

    ! create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text

    ! allocate scalars
    call pcpobj%pcp_allocate_scalars()

    ! initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'DFW')

    ! store pointer to dis
    pcpobj%dis => dis

    ! store pointer to dfw
    pcpobj%dfw => dfw

    ! store pointer to cxs
    pcpobj%cxs => cxs
  end subroutine pcp_create

  !> @brief Allocate scalar members
  !<
  subroutine pcp_allocate_scalars(this)
    ! dummy
    class(SwfPcpType), intent(inout) :: this

    ! allocate base scalars
    call this%BndExtType%allocate_scalars()

    ! allocate internal members
    allocate (this%read_as_arrays)

    ! Set values
    this%read_as_arrays = .false.
  end subroutine pcp_allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine pcp_allocate_arrays(this, nodelist, auxvar)
    ! modules
    use MemoryManagerModule, only: mem_setptr, mem_checkin
    ! dummy
    class(SwfPcpType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar

    ! allocate base arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)

    ! set input context pointers
    call mem_setptr(this%precipitation, 'PRECIPITATION', this%input_mempath)

    ! checkin input context pointers
    call mem_checkin(this%precipitation, 'PRECIPITATION', this%memoryPath, &
                     'PRECIPITATION', this%input_mempath)
  end subroutine pcp_allocate_arrays

  !> @brief Source options specific to PCPType
  !<
  subroutine pcp_source_options(this)
    ! modules
    use MemoryManagerExtModule, only: mem_set_value
    implicit none
    ! dummy
    class(SwfPcpType), intent(inout) :: this
    ! local
    logical(LGP) :: found_readasarrays = .false.

    ! source common bound options
    call this%BndExtType%source_options()

    ! update defaults with idm sourced values
    call mem_set_value(this%read_as_arrays, 'READASARRAYS', this%input_mempath, &
                       found_readasarrays)

    ! log pcp params
    call this%log_pcp_options(found_readasarrays)
  end subroutine pcp_source_options

  !> @brief Log options specific to SwfPcpType
  !<
  subroutine log_pcp_options(this, found_readasarrays)
    implicit none
    ! dummy
    class(SwfPcpType), intent(inout) :: this
    logical(LGP), intent(in) :: found_readasarrays
    ! formats
    character(len=*), parameter :: fmtreadasarrays = &
      &"(4x, 'PRECIPITATION INPUT WILL BE READ AS ARRAY(S).')"

    ! log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'

    if (found_readasarrays) then
      write (this%iout, fmtreadasarrays)
    end if

    ! close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine log_pcp_options

  !> @brief Source the dimensions for this package
  !<
  subroutine pcp_source_dimensions(this)
    ! dummy
    class(SwfPcpType), intent(inout) :: this

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
  end subroutine pcp_source_dimensions

  !> @brief Part of allocate and read
  !<
  subroutine pcp_read_initial_attr(this)
    ! dummy
    class(SwfPcpType), intent(inout) :: this

    if (this%read_as_arrays) then
      call this%default_nodelist()
    end if
  end subroutine pcp_read_initial_attr

  !> @brief Read and Prepare
  !!
  !! Read itmp and read new boundaries if itmp > 0
  !<
  subroutine pcp_rp(this)
    ! modules
    use TdisModule, only: kper
    implicit none
    ! dummy
    class(SwfPcpType), intent(inout) :: this

    if (this%iper /= kper) return

    if (this%read_as_arrays) then
      ! no need to do anything because this%precipitation points directly to
      ! the input context precipitation, which is automatically updated by idm
    else
      call this%BndExtType%bnd_rp()
    end if

    ! Write the list to iout if requested
    if (this%iprpak /= 0) then
      call this%write_list()
    end if
  end subroutine pcp_rp

  !> @brief Ensure precipitation is positive
  !<
  subroutine pcp_ck(this)
    ! dummy
    class(SwfPcpType), intent(inout) :: this
    ! local
    character(len=30) :: nodestr
    integer(I4B) :: i, nr
    character(len=*), parameter :: fmterr = &
      &"('Specified stress ',i0, &
      &' precipitation (',g0,') is less than zero for cell', a)"

    ! Ensure precipitation rates are positive
    do i = 1, this%nbound
      nr = this%nodelist(i)
      if (nr <= 0) cycle
      if (this%precipitation(i) < DZERO) then
        call this%dis%noder_to_string(nr, nodestr)
        write (errmsg, fmt=fmterr) i, this%precipitation(i), trim(nodestr)
        call store_error(errmsg)
      end if
    end do

    ! write summary of package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine pcp_ck

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Skip if no precipitation. Otherwise, calculate hcof and rhs
  !<
  subroutine pcp_cf(this)
    ! dummy
    class(SwfPcpType) :: this
    ! local
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: idcxs
    real(DP) :: qpcp
    real(DP) :: area
    real(DP) :: width_channel
    real(DP) :: top_width
    real(DP) :: dummy
    real(DP), dimension(:), pointer :: reach_length

    ! Return if no precipitation
    if (this%nbound == 0) return

    ! Set pointer to reach_length for 1d
    reach_length => this%reach_length_pointer()

    ! Calculate hcof and rhs for each precipitation entry
    do i = 1, this%nbound

      ! Find the node number
      node = this%nodelist(i)

      ! cycle if nonexistent bound
      if (node <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      ! Initialize hcof
      this%hcof(i) = DZERO

      ! Determine the water surface area
      if (this%dis%is_2d()) then
        ! this is for overland flow case
        area = this%dis%get_area(node)
      else if (this%dis%is_1d()) then
        ! this is for channel case
        idcxs = this%dfw%idcxs(node)
        call this%dis%get_flow_width(node, node, 0, width_channel, &
                                     dummy)
        top_width = this%cxs%get_maximum_top_width(idcxs, width_channel)
        area = reach_length(node) * top_width
      end if

      ! calculate volumetric precipitation flow in L^3/T
      qpcp = this%precipitation(i) * area

      ! multiplier
      if (this%iauxmultcol > 0) then
        qpcp = qpcp * this%auxvar(this%iauxmultcol, i)
      end if

      ! rhs contribution
      this%rhs(i) = -qpcp

      ! zero out contribution if cell is inactive or constant head
      if (this%ibound(node) <= 0) then
        this%rhs(i) = DZERO
        cycle
      end if

    end do
  end subroutine pcp_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine pcp_fc(this, rhs, ia, idxglo, matrix_sln)
    ! dummy
    class(SwfPcpType) :: this
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
  end subroutine pcp_fc

  !> @brief Deallocate memory
  !<
  subroutine pcp_da(this)
    ! modules
    use MemoryManagerModule, only: mem_deallocate
    ! dummy
    class(SwfPcpType) :: this

    ! Deallocate parent package
    call this%BndExtType%bnd_da()

    ! scalars
    deallocate (this%read_as_arrays)

    ! arrays
    call mem_deallocate(this%precipitation, 'PRECIPITATION', this%memoryPath)

    ! pointers
    nullify (this%dis)
    nullify (this%dfw)
    nullify (this%cxs)
  end subroutine pcp_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used.
  !<
  subroutine pcp_define_listlabel(this)
    ! dummy
    class(SwfPcpType), intent(inout) :: this
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'PRECIPITATION'
!    if(this%multindex > 0) &
!      write(this%listlabel, '(a, a16)') trim(this%listlabel), 'MULTIPLIER'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine pcp_define_listlabel

  !> @brief Assign default nodelist when READASARRAYS is specified.
  !<
  subroutine default_nodelist(this)
    ! dummy
    class(SwfPcpType) :: this
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
  logical function pcp_obs_supported(this)
    implicit none
    ! dummy
    class(SwfPcpType) :: this
    pcp_obs_supported = .true.
  end function pcp_obs_supported

  !> @brief Implements bnd_df_obs
  !!
  !! Store observation type supported by PCP package. Overrides
  !! BndType%bnd_df_obs
  !<
  subroutine pcp_df_obs(this)
    implicit none
    ! dummy
    class(SwfPcpType) :: this
    ! local
    integer(I4B) :: indx

    call this%obs%StoreObsType('pcp', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine pcp_df_obs

  !> @brief Return requested boundary value
  !<
  function pcp_bound_value(this, col, row) result(bndval)
    ! modules
    use ConstantsModule, only: DZERO
    ! dummy
    class(SwfPcpType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! result
    real(DP) :: bndval

    select case (col)
    case (1)
      if (this%iauxmultcol > 0) then
        bndval = this%precipitation(row) * this%auxvar(this%iauxmultcol, row)
      else
        bndval = this%precipitation(row)
      end if
    case default
      errmsg = 'Programming error. PCP bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function pcp_bound_value

  function reach_length_pointer(this) result(ptr)
    ! dummy
    class(SwfPcpType) :: this !< this instance
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

end module SwfPcpModule

