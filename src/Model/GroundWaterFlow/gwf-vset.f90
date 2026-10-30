module GwfVsetModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DHALF, DONE, LENFTYPE, LENPACKAGENAME, &
                             LENVARNAME, LENMEMPATH, IWETLAKE
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use SimModule, only: store_error, store_error_filename, count_errors
  use SimVariablesModule, only: errmsg
  use ObsModule, only: DefaultObsIdProcessor
  use MatrixBaseModule
  use EtUptakeModelModule, only: EtUptakeModelType
  use EtUptakeModelFactoryModule, only: create_uptake_model, &
                                        destroy_uptake_model, &
                                        uptake_model_name, FEDDES
  implicit none
  private

  public :: vset_create

  character(len=LENFTYPE) :: ftype = 'VSET'
  character(len=LENPACKAGENAME) :: text = '            VSET'

  type, extends(BndExtType) :: GwfVsetType
    integer(I4B), pointer :: et_model_id => null() !< selected uptake model id
    real(DP), pointer :: h1 => null() !< Feddes anoxia pressure head
    real(DP), pointer :: h2 => null() !< Feddes wet end of optimal range
    real(DP), pointer :: h3 => null() !< Feddes dry end of optimal range
    real(DP), pointer :: h4 => null() !< Feddes wilting point pressure head
    real(DP), dimension(:), pointer, contiguous :: rate => null() !< potential transpiration rate
    integer(I4B), dimension(:), pointer, contiguous :: iunsat => null() !< VSF unsaturated-cell indicator (not owned)
    class(EtUptakeModelType), pointer :: uptake_model => null() !< pressure-head uptake model
  contains
    procedure :: allocate_scalars => vset_allocate_scalars
    procedure :: allocate_arrays => vset_allocate_arrays
    procedure :: source_options => vset_source_options
    procedure :: bnd_rp => vset_rp
    procedure :: bnd_cf => vset_cf
    procedure :: bnd_fc => vset_fc
    procedure :: bnd_da => vset_da
    procedure :: define_listlabel => vset_define_listlabel
    procedure :: bound_value => vset_bound_value
    procedure, public :: bnd_obs_supported => vset_obs_supported
    procedure, public :: bnd_df_obs => vset_df_obs
    ! private
    procedure, private :: vset_log_options
    procedure, private :: setup_uptake_model => vset_setup_uptake_model
  end type GwfVsetType

contains

  !> @brief Create a new Variably-Saturated Evapotranspiration Package
  !<
  subroutine vset_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                         mempath)
    class(BndType), pointer :: packobj !< the new package
    integer(I4B), intent(in) :: id !< consecutive package number
    integer(I4B), intent(in) :: ibcnum !< consecutive boundary package number
    integer(I4B), intent(in) :: inunit !< the input file unit
    integer(I4B), intent(in) :: iout !< the output file unit
    character(len=*), intent(in) :: namemodel !< the model name
    character(len=*), intent(in) :: pakname !< the package name
    character(len=*), intent(in) :: mempath !< the input memory path
    ! local
    type(GwfVsetType), pointer :: vsetobj

    allocate (vsetobj)
    packobj => vsetobj

    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text

    call vsetobj%allocate_scalars()
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')

  end subroutine vset_create

  !> @brief Allocate package scalars
  !<
  subroutine vset_allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    class(GwfVsetType) :: this !< this instance

    call this%BndExtType%allocate_scalars()

    call mem_allocate(this%et_model_id, 'ET_MODEL_ID', this%memoryPath)
    call mem_allocate(this%h1, 'H1', this%memoryPath)
    call mem_allocate(this%h2, 'H2', this%memoryPath)
    call mem_allocate(this%h3, 'H3', this%memoryPath)
    call mem_allocate(this%h4, 'H4', this%memoryPath)

    this%et_model_id = FEDDES
    this%h1 = DZERO
    this%h2 = DZERO
    this%h3 = DZERO
    this%h4 = DZERO

  end subroutine vset_allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine vset_allocate_arrays(this, nodelist, auxvar)
    use MemoryManagerModule, only: mem_setptr, mem_checkin, get_isize
    class(GwfVsetType) :: this !< this instance
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist !< package nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar !< package aux variables
    ! local
    character(len=LENMEMPATH) :: vsf_mempath
    integer(I4B) :: isize

    call this%BndExtType%allocate_arrays(nodelist, auxvar)

    ! set and checkin input context pointers
    call mem_setptr(this%rate, 'RATE', this%input_mempath)
    call mem_checkin(this%rate, 'RATE', this%memoryPath, &
                     'RATE', this%input_mempath)

    ! VSET can only be applied on VSF (Richards) cells, so VSF must be active;
    ! bind its unsaturated-cell indicator for per-cell validation
    vsf_mempath = create_mem_path(this%name_model, 'VSF')
    call get_isize('IUNSAT', vsf_mempath, isize)
    if (isize < 0) then
      write (errmsg, '(a)') 'The VSET package requires the VSF package to be '// &
        'active in the same GWF model.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
      return
    end if
    call mem_setptr(this%iunsat, 'IUNSAT', vsf_mempath)

  end subroutine vset_allocate_arrays

  !> @brief Source options specific to the VSET package
  !<
  subroutine vset_source_options(this)
    use MemoryManagerExtModule, only: mem_set_value
    use GwfVsetInputModule, only: GwfVsetParamFoundType
    class(GwfVsetType), intent(inout) :: this !< this instance
    ! local
    type(GwfVsetParamFoundType) :: found
    character(len=LENVARNAME) :: et_model_str

    call this%BndExtType%source_options()

    et_model_str = ''
    call mem_set_value(et_model_str, 'ET_MODEL', this%input_mempath, &
                       found%et_model)
    call mem_set_value(this%h1, 'H1', this%input_mempath, found%h1)
    call mem_set_value(this%h2, 'H2', this%input_mempath, found%h2)
    call mem_set_value(this%h3, 'H3', this%input_mempath, found%h3)
    call mem_set_value(this%h4, 'H4', this%input_mempath, found%h4)

    ! resolve the uptake model
    if (found%et_model) then
      select case (trim(et_model_str))
      case ('FEDDES')
        this%et_model_id = FEDDES
      case default
        write (errmsg, '(a,a)') 'Unknown ET_MODEL: ', trim(et_model_str)
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end select
    end if

    call this%vset_log_options(found)
    call this%setup_uptake_model()

  end subroutine vset_source_options

  !> @brief Log VSET specific options
  !<
  subroutine vset_log_options(this, found)
    use GwfVsetInputModule, only: GwfVsetParamFoundType
    class(GwfVsetType), intent(inout) :: this !< this instance
    type(GwfVsetParamFoundType), intent(in) :: found !< sourced input flags

    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'

    write (this%iout, '(4x,a,a)') 'ET MODEL IS ', &
      trim(uptake_model_name(this%et_model_id))

    if (found%h1) write (this%iout, '(4x,a,g15.7)') 'FEDDES H1 = ', this%h1
    if (found%h2) write (this%iout, '(4x,a,g15.7)') 'FEDDES H2 = ', this%h2
    if (found%h3) write (this%iout, '(4x,a,g15.7)') 'FEDDES H3 = ', this%h3
    if (found%h4) write (this%iout, '(4x,a,g15.7)') 'FEDDES H4 = ', this%h4

    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'

  end subroutine vset_log_options

  !> @brief Validate parameters and create the uptake model instance
  !<
  subroutine vset_setup_uptake_model(this)
    use EtUptakeFeddesModule, only: FeddesUptakeType
    class(GwfVsetType), intent(inout) :: this !< this instance

    select case (this%et_model_id)
    case (FEDDES)
      ! thresholds must be monotonically decreasing h1 > h2 > h3 > h4
      if (.not. (this%h1 > this%h2 .and. this%h2 > this%h3 .and. &
                 this%h3 > this%h4)) then
        write (errmsg, '(a)') 'FEDDES pressure-head thresholds must satisfy '// &
          'H1 > H2 > H3 > H4.  Specify all four in the OPTIONS block.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
    end select

    if (count_errors() > 0) return

    this%uptake_model => create_uptake_model(this%et_model_id)
    select type (um => this%uptake_model)
    type is (FeddesUptakeType)
      um%h1 = this%h1
      um%h2 = this%h2
      um%h3 = this%h3
      um%h4 = this%h4
    end select

  end subroutine vset_setup_uptake_model

  !> @brief Read and prepare; validate that each cell is a VSF cell
  !<
  subroutine vset_rp(this)
    class(GwfVsetType), intent(inout) :: this !< this instance
    ! local
    integer(I4B) :: i, node

    call this%BndExtType%bnd_rp()

    if (.not. associated(this%iunsat)) return

    do i = 1, this%nbound
      node = this%nodelist(i)
      if (node <= 0) cycle
      if (this%iunsat(node) /= 1) then
        write (errmsg, '(a,i0,a)') 'VSET cell at reduced node ', node, &
          ' is not a VSF (Richards) cell.  VSET may only be applied where '// &
          'IUNSAT is active in the VSF package.'
        call store_error(errmsg)
      end if
    end do
    if (count_errors() > 0) call store_error_filename(this%input_fname)

  end subroutine vset_rp

  !> @brief Formulate the HCOF and RHS terms
  !<
  subroutine vset_cf(this)
    class(GwfVsetType) :: this !< this instance
    ! local
    integer(I4B) :: i, node
    real(DP) :: c, h, z, hcof, rhs

    if (this%nbound == 0) return

    do i = 1, this%nbound
      node = this%nodelist(i)

      this%rhs(i) = DZERO
      this%hcof(i) = DZERO

      if (node <= 0) cycle

      ! skip inactive or lake-overlain cells
      if (this%ibound(node) <= 0 .or. this%ibound(node) == IWETLAKE) cycle

      ! potential (volumetric) transpiration rate for the cell
      if (this%iauxmultcol > 0) then
        c = this%rate(i) * this%dis%get_area(node) * &
            this%auxvar(this%iauxmultcol, i)
      else
        c = this%rate(i) * this%dis%get_area(node)
      end if

      ! cell-center elevation and current head
      z = DHALF * (this%dis%bot(node) + this%dis%top(node))
      h = this%xnew(node)

      ! linearized uptake terms from the selected model
      call this%uptake_model%fc_terms(h, z, c, hcof, rhs)
      this%hcof(i) = hcof
      this%rhs(i) = rhs
    end do

  end subroutine vset_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine vset_fc(this, rhs, ia, idxglo, matrix_sln)
    class(GwfVsetType) :: this !< this instance
    real(DP), dimension(:), intent(inout) :: rhs !< solution right-hand side
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< global matrix indices
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! local
    integer(I4B) :: i, n, ipos

    do i = 1, this%nbound
      n = this%nodelist(i)
      if (n <= 0) cycle
      ! reset hcof and rhs for excluded cells
      if (this%ibound(n) == IWETLAKE) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
    end do

  end subroutine vset_fc

  !> @brief Clean up
  !<
  subroutine vset_da(this)
    use MemoryManagerModule, only: mem_deallocate
    class(GwfVsetType) :: this !< this instance

    call destroy_uptake_model(this%uptake_model)

    ! iunsat is owned by the VSF package; just drop the alias
    this%iunsat => null()

    call mem_deallocate(this%et_model_id)
    call mem_deallocate(this%h1)
    call mem_deallocate(this%h2)
    call mem_deallocate(this%h3)
    call mem_deallocate(this%h4)

    call mem_deallocate(this%rate, 'RATE', this%memoryPath)

    call this%BndExtType%bnd_da()

  end subroutine vset_da

  !> @brief Define the list heading written to iout when PRINT_INPUT is used
  !<
  subroutine vset_define_listlabel(this)
    class(GwfVsetType), intent(inout) :: this !< this instance

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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'POT. RATE'

    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if

  end subroutine vset_define_listlabel

  !> @brief Return requested boundary value
  !<
  function vset_bound_value(this, col, row) result(bndval)
    class(GwfVsetType), intent(inout) :: this !< this instance
    integer(I4B), intent(in) :: col !< bound column
    integer(I4B), intent(in) :: row !< bound row
    real(DP) :: bndval !< the returned bound value

    select case (col)
    case (1)
      if (this%iauxmultcol > 0) then
        bndval = this%rate(row) * this%auxvar(this%iauxmultcol, row)
      else
        bndval = this%rate(row)
      end if
    case default
      write (errmsg, '(a,i0,a)') &
        'Programming error. VSET bound value requested for column ', col, '.'
      call store_error(errmsg, terminate=.TRUE.)
    end select

  end function vset_bound_value

  !> @brief Return true because the VSET package supports observations
  !<
  logical function vset_obs_supported(this)
    class(GwfVsetType) :: this !< this instance

    vset_obs_supported = .true.

  end function vset_obs_supported

  !> @brief Store observation type supported by the VSET package
  !<
  subroutine vset_df_obs(this)
    class(GwfVsetType) :: this !< this instance
    ! local
    integer(I4B) :: indx

    call this%obs%StoreObsType('vset', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor

  end subroutine vset_df_obs

end module GwfVsetModule
