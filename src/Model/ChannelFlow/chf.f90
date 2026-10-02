!> @brief Channel Flow (CHF) Module
!<
module ChfModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LENPACKAGETYPE, LENMEMPATH, LINELENGTH
  use SimModule, only: store_error
  use BaseModelModule, only: BaseModelType
  use ListsModule, only: basemodellist
  use BaseModelModule, only: AddBaseModelToList
  use BndModule, only: BndType, GetBndFromList
  use SwfModule, only: SwfModelType, swf_ac, swf_mc, swf_fc, swf_cq, swf_da
  use BudgetModule, only: budget_cr
  use MatrixBaseModule, only: MatrixBaseType
  use SparseModule, only: sparsematrix
  use ChfJncModule, only: ChfJncType, chf_jnc_cr

  implicit none

  private
  public :: chf_cr
  public :: ChfModelType
  public :: CHF_NBASEPKG, CHF_NMULTIPKG
  public :: CHF_BASEPKG, CHF_MULTIPKG

  type, extends(SwfModelType) :: ChfModelType
    type(ChfJncType), pointer :: jnc => null() !< channel junction package (CHF-only)
  contains
    procedure :: model_df => chf_df
    procedure :: model_ac => chf_ac
    procedure :: model_mc => chf_mc
    procedure :: model_fc => chf_fc
    procedure :: model_cq => chf_cq
    procedure :: model_da => chf_da
    procedure :: set_namfile_options
    procedure :: log_namfile_options
  end type ChfModelType

  !> @brief CHF base package array descriptors
  !!
  !! CHF model base package types.  Only listed packages are candidates
  !< for input and these will be loaded in the order specified.
  integer(I4B), parameter :: CHF_NBASEPKG = 7
  character(len=LENPACKAGETYPE), dimension(CHF_NBASEPKG) :: &
    CHF_BASEPKG = ['DISV1D6', 'DFW6   ', 'CXS6   ', &
                   'OC6    ', 'IC6    ', 'OBS6   ', &
                   'STO6   ']

  !> @brief CHF multi package array descriptors
  !!
  !! CHF model multi-instance package types.  Only listed packages are
  !< candidates for input and these will be loaded in the order specified.
  integer(I4B), parameter :: CHF_NMULTIPKG = 50
  character(len=LENPACKAGETYPE), dimension(CHF_NMULTIPKG) :: CHF_MULTIPKG
  data CHF_MULTIPKG/'FLW6 ', 'CHD6 ', 'CDB6 ', 'ZDG6 ', 'PCP6 ', & !  5
                    'EVP6 ', '     ', '     ', '     ', '     ', & !  10
                    &40*'     '/ ! 50

  ! size of supported model package arrays
  integer(I4B), parameter :: NIUNIT_CHF = CHF_NBASEPKG + CHF_NMULTIPKG

contains

  !> @brief Create a new surface water flow model object
  !<
  subroutine chf_cr(filename, id, modelname)
    ! modules
    ! dummy
    character(len=*), intent(in) :: filename !< input file
    integer(I4B), intent(in) :: id !< consecutive model number listed in mfsim.nam
    character(len=*), intent(in) :: modelname !< name of the model
    ! local
    class(ChfModelType), pointer :: this
    class(BaseModelType), pointer :: model

    ! Allocate a new model
    allocate (this)
    model => this
    call AddBaseModelToList(basemodellist, model)

    ! call parent initialize routine
    call this%initialize('CHF', filename, id, modelname)

    ! set and log namefile options
    call this%set_namfile_options()

    ! Create utility objects
    call budget_cr(this%budget, this%name)

    ! create model packages
    call this%create_packages()

    ! create the CHF junction package (auto-generated, CHF-only)
    call chf_jnc_cr(this%jnc, this%name, this%iout)

  end subroutine chf_cr

  !> @brief Define the CHF model
  !!
  !! Overrides SwfModelType%model_df (swf_df) to detect channel junctions and
  !! grow the model equation count to include them, before model arrays are
  !! allocated.  The junction rows are appended after the reach rows; their
  !! connections are added to the solution sparse matrix in chf_ac (the model's
  !! own ia/ja, which point at dis%con, are left untouched - MAW convention).
  !<
  subroutine chf_df(this)
    ! dummy
    class(ChfModelType) :: this
    ! local
    integer(I4B) :: ip
    class(BndType), pointer :: packobj

    ! call package df routines
    call this%dis%dis_df()
    call this%dfw%dfw_df(this%dis)
    call this%oc%oc_df()
    call this%budget%budget_df(NIUNIT_CHF, 'VOLUME', 'L**3')

    ! detect junctions from the DISV1D connectivity and mask the direct
    ! reach-reach connections they replace (so DFW routes through the junction)
    call this%jnc%set_pointers(this%dis, this%dfw)
    call this%jnc%detect_junctions()
    call this%jnc%mask_reach_connections()

    ! junction rows are appended immediately after the reach rows
    this%jnc%ioffset = 0

    ! set model sizes: reach equations plus one per junction.  nja and ia/ja
    ! stay as the DIS reach structure; junction connections are added only to
    ! the solution sparse matrix (chf_ac), never to dis%con.
    this%neq = this%dis%nodes + this%jnc%njunctions
    this%nja = this%dis%nja
    this%ia => this%dis%con%ia
    this%ja => this%dis%con%ja

    ! allocate model arrays, now that neq and nja are known
    call this%allocate_arrays()

    ! define packages and assign iout for time series managers
    do ip = 1, this%bndlist%Count()
      packobj => GetBndFromList(this%bndlist, ip)
      call packobj%bnd_df(this%dis%nodes, this%dis)
    end do

    ! store information needed for observations
    call this%obs%obs_df(this%iout, this%name, 'SWF', this%dis)

  end subroutine chf_df

  !> @brief Add CHF model connections to the sparse matrix
  !!
  !! Adds the reach grid connections (shared swf_ac) plus the junction rows and
  !! reach-junction connections (jnc_ac).
  !<
  subroutine chf_ac(this, sparse)
    ! dummy
    class(ChfModelType) :: this
    type(sparsematrix), intent(inout) :: sparse

    ! reach connections and any boundary-package connections
    call swf_ac(this, sparse)

    ! junction rows and reach-junction connections
    call this%jnc%jnc_ac(this%moffset, this%dis%nodes, sparse)

  end subroutine chf_ac

  !> @brief Map CHF model connection positions in the solution matrix
  !<
  subroutine chf_mc(this, matrix_sln)
    ! dummy
    class(ChfModelType) :: this
    class(MatrixBaseType), pointer :: matrix_sln

    ! reach and boundary-package position mapping
    call swf_mc(this, matrix_sln)

    ! junction connection position caching
    call this%jnc%jnc_mc(this%moffset, this%dis%nodes, matrix_sln)

  end subroutine chf_mc

  !> @brief Fill CHF model coefficients
  !!
  !! Runs the shared SWF fill (DFW, storage, boundary packages) then adds the
  !! junction continuity rows.
  !<
  subroutine chf_fc(this, kiter, matrix_sln, inwtflag)
    ! dummy
    class(ChfModelType) :: this
    integer(I4B), intent(in) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in) :: inwtflag

    ! shared SWF coefficient fill
    call swf_fc(this, kiter, matrix_sln, inwtflag)

    ! junction continuity rows
    call this%jnc%jnc_fc(matrix_sln, this%rhs, this%x, this%dis%nodes)

  end subroutine chf_fc

  !> @brief Calculate CHF model flows
  !<
  subroutine chf_cq(this, icnvg, isuppress_output)
    ! dummy
    class(ChfModelType) :: this
    integer(I4B), intent(in) :: icnvg
    integer(I4B), intent(in) :: isuppress_output

    ! shared SWF flow calculation
    call swf_cq(this, icnvg, isuppress_output)

    ! junction stage capture and reach-junction flow contribution to flowja
    call this%jnc%jnc_cq(this%x, this%flowja, this%dis%nodes)

  end subroutine chf_cq

  !> @brief Deallocate the CHF model
  !<
  subroutine chf_da(this)
    ! dummy
    class(ChfModelType) :: this

    ! deallocate the junction package
    call this%jnc%jnc_da()
    deallocate (this%jnc)
    nullify (this%jnc)

    ! deallocate the shared SWF model
    call swf_da(this)

  end subroutine chf_da

  !> @brief Handle namefile options
  !!
  !! Set pointers to IDM namefile options, then
  !! create the list file and log options.
  !<
  subroutine set_namfile_options(this)
    use SimVariablesModule, only: idm_context
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerExtModule, only: mem_set_value
    use ChfNamInputModule, only: ChfNamParamFoundType
    class(ChfModelType) :: this
    type(ChfNamParamFoundType) :: found
    character(len=LENMEMPATH) :: input_mempath
    character(len=LINELENGTH) :: lst_fname

    ! set input model namfile memory path
    input_mempath = create_mem_path(this%name, 'NAM', idm_context)

    ! copy option params from input context
    call mem_set_value(lst_fname, 'LIST', input_mempath, found%list)
    call mem_set_value(this%inewton, 'NEWTON', input_mempath, found%newton)
    call mem_set_value(this%inewtonur, 'UNDER_RELAXATION', input_mempath, &
                       found%under_relaxation)
    call mem_set_value(this%iprpak, 'PRINT_INPUT', input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', input_mempath, &
                       found%print_flows)
    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', input_mempath, found%save_flows)

    ! create the list file
    call this%create_lstfile(lst_fname, this%filename, found%list, &
                             'CHANNEL FLOW MODEL (CHF)')

    ! activate save_flows if found
    if (found%save_flows) then
      this%ipakcb = -1
    end if

    ! log set options
    if (this%iout > 0) then
      call this%log_namfile_options(found)
    end if

  end subroutine set_namfile_options

  !> @brief Write model namfile options to list file
  !<
  subroutine log_namfile_options(this, found)
    use ChfNamInputModule, only: ChfNamParamFoundType
    class(ChfModelType) :: this
    type(ChfNamParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'BEGIN NAMEFILE OPTIONS'

    if (found%print_input) then
      write (this%iout, '(4x,a)') 'STRESS PACKAGE INPUT WILL BE PRINTED '// &
        'FOR ALL MODEL STRESS PACKAGES'
    end if

    if (found%print_flows) then
      write (this%iout, '(4x,a)') 'PACKAGE FLOWS WILL BE PRINTED '// &
        'FOR ALL MODEL PACKAGES'
    end if

    if (found%save_flows) then
      write (this%iout, '(4x,a)') &
        'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL'
    end if

    write (this%iout, '(1x,a)') 'END NAMEFILE OPTIONS'

  end subroutine log_namfile_options

end module ChfModule
