module GwtDspModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DHALF, DPI
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use Xt3dModule, only: Xt3dType, xt3d_cr
  use GwtDspOptionsModule, only: GwtDspOptionsType
  use MatrixBaseModule

  implicit none
  private
  public :: GwtDspType
  public :: dsp_cr

  type, extends(NumericalPackageType) :: GwtDspType

    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() ! pointer to GWT model ibound
    type(TspFmiType), pointer :: fmi => null() ! pointer to GWT fmi object
    real(DP), dimension(:), pointer, contiguous :: thetam => null() ! pointer to GWT storage porosity (voids per aquifer volume)
    real(DP), dimension(:), pointer, contiguous :: diffc => null() ! molecular diffusion coefficient for each cell
    real(DP), dimension(:), pointer, contiguous :: alh => null() ! longitudinal horizontal dispersivity
    real(DP), dimension(:), pointer, contiguous :: alv => null() ! longitudinal vertical dispersivity
    real(DP), dimension(:), pointer, contiguous :: ath1 => null() ! transverse horizontal dispersivity
    real(DP), dimension(:), pointer, contiguous :: ath2 => null() ! transverse horizontal dispersivity
    real(DP), dimension(:), pointer, contiguous :: atv => null() ! transverse vertical dispersivity
    integer(I4B), pointer :: idiffc => null() ! flag indicating diffusion is active
    integer(I4B), pointer :: idisp => null() ! flag indicating mechanical dispersion is active
    integer(I4B), pointer :: ialh => null() ! longitudinal horizontal dispersivity data flag
    integer(I4B), pointer :: ialv => null() ! longitudinal vertical dispersivity data flag
    integer(I4B), pointer :: iath1 => null() ! transverse horizontal dispersivity data flag
    integer(I4B), pointer :: iath2 => null() ! transverse horizontal dispersivity data flag
    integer(I4B), pointer :: iatv => null() ! transverse vertical dispersivity data flag
    integer(I4B), pointer :: ixt3doff => null() ! xt3d off flag, xt3d is set inactive if 1
    integer(I4B), pointer :: ixt3drhs => null() ! xt3d rhs flag, xt3d rhs is set active if 1
    integer(I4B), pointer :: ixt3d => null() ! flag indicating xt3d is active
    type(Xt3dType), pointer :: xt3d => null() ! xt3d object
    real(DP), dimension(:), pointer, contiguous :: dispcoef => null() ! disp coefficient (only if xt3d not active)
    integer(I4B), pointer :: id22 => null() ! flag indicating d22 is available
    integer(I4B), pointer :: id33 => null() ! flag indicating d33 is available
    real(DP), dimension(:), pointer, contiguous :: d11 => null() ! dispersion coefficient
    real(DP), dimension(:), pointer, contiguous :: d22 => null() ! dispersion coefficient
    real(DP), dimension(:), pointer, contiguous :: d33 => null() ! dispersion coefficient
    real(DP), dimension(:), pointer, contiguous :: angle1 => null() ! rotation angle 1
    real(DP), dimension(:), pointer, contiguous :: angle2 => null() ! rotation angle 2
    real(DP), dimension(:), pointer, contiguous :: angle3 => null() ! rotation angle 3
    integer(I4B), pointer :: iangle1 => null() ! flag indicating angle1 is available
    integer(I4B), pointer :: iangle2 => null() ! flag indicating angle2 is available
    integer(I4B), pointer :: iangle3 => null() ! flag indicating angle3 is available

  contains

    procedure :: dsp_df
    procedure :: dsp_ac
    procedure :: dsp_mc
    procedure :: dsp_ar
    procedure :: dsp_ad
    procedure :: dsp_fc
    procedure :: dsp_cq
    procedure :: dsp_da
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_griddata
    procedure, private :: log_options
    procedure, private :: log_griddata
    procedure, private :: calcdispellipse
    procedure, private :: calcdispcoef

  end type GwtDspType

contains

  !> @brief Create a DSP object
  !<
  subroutine dsp_cr(dspobj, name_model, input_mempath, inunit, iout, fmi)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    type(GwtDspType), pointer :: dspobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi
    ! -- locals
    ! -- formats
    character(len=*), parameter :: fmtdsp = &
      "(1x,/1x,'DSP-- DISPERSION PACKAGE, VERSION 1, 1/24/2018', &
      &' INPUT READ FROM MEMPATH: ', A, //)"
    !
    ! -- Create the object
    allocate (dspobj)
    !
    ! -- create name and memory path
    call dspobj%set_names(1, name_model, 'DSP', 'DSP', input_mempath)
    !
    ! -- Allocate scalars
    call dspobj%allocate_scalars()
    !
    ! -- Set variables
    dspobj%inunit = inunit
    dspobj%iout = iout
    dspobj%fmi => fmi
    !
    if (dspobj%inunit > 0) then
      !
      ! -- Print a message identifying the dispersion package.
      if (dspobj%iout > 0) then
        write (dspobj%iout, fmtdsp) input_mempath
      end if
    end if
  end subroutine dsp_cr

  !> @brief Define DSP object
  !<
  subroutine dsp_df(this, dis, dspOptions)
    ! -- modules
    ! -- dummy
    class(GwtDspType) :: this
    class(DisBaseType), pointer :: dis
    type(GwtDspOptionsType), optional, intent(in) :: dspOptions !< the optional DSP options, used when not
                                                                !! creating DSP from file
    ! -- local
    !
    ! -- Store pointer to dis
    this%dis => dis
    !
    !
    ! -- set default xt3d representation to on and lhs
    this%ixt3d = 1
    !
    ! -- Read dispersion options
    if (present(dspOptions)) then
      this%ixt3d = dspOptions%ixt3d
      !
      ! -- Allocate only, grid data will not be read from file
      call this%allocate_arrays(this%dis%nodes)
    else
      !
      ! -- Source options
      call this%source_options()
      call this%allocate_arrays(this%dis%nodes)
      !
      ! -- Source dispersion data
      call this%source_griddata()
    end if
    !
    ! -- xt3d create
    if (this%ixt3d > 0) then
      call xt3d_cr(this%xt3d, this%name_model, this%inunit, this%iout, &
                   ldispopt=.true.)
      this%xt3d%ixt3d = this%ixt3d
      call this%xt3d%xt3d_df(dis)
    end if
  end subroutine dsp_df

  !> @brief Add connections to DSP
  !!
  !! Add connections for extended neighbors to the sparse matrix
  !<
  subroutine dsp_ac(this, moffset, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtDspType) :: this
    integer(I4B), intent(in) :: moffset
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    !
    ! -- Add extended neighbors (neighbors of neighbors)
    if (this%ixt3d > 0) call this%xt3d%xt3d_ac(moffset, sparse)
  end subroutine dsp_ac

  !> @brief Map DSP connections
  !!
  !! Map connections and construct iax, jax, and idxglox
  !<
  subroutine dsp_mc(this, moffset, matrix_sln)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtDspType) :: this
    integer(I4B), intent(in) :: moffset
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    !
    ! -- Call xt3d map connections
    if (this%ixt3d > 0) call this%xt3d%xt3d_mc(moffset, matrix_sln)
  end subroutine dsp_mc

  !> @brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the package.
  !<
  subroutine dsp_ar(this, ibound, thetam)
    ! -- modules
    ! -- dummy
    class(GwtDspType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    real(DP), dimension(:), pointer, contiguous :: thetam
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtdsp = &
      "(1x,/1x,'DSP-- DISPERSION PACKAGE, VERSION 1, 1/24/2018', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- dsp pointers to arguments that were passed in
    this%ibound => ibound
    this%thetam => thetam
  end subroutine dsp_ar

  !> @brief Advance method for the package
  !<
  subroutine dsp_ad(this)
    ! -- modules
    use TdisModule, only: kstp, kper
    ! -- dummy
    class(GwtDspType) :: this
    ! -- local
    !
    ! -- xt3d
    ! TODO: might consider adding a new mf6 level set pointers method, and
    ! doing this stuff there instead of in the time step loop.
    if (kstp * kper == 1) then
      if (this%ixt3d > 0) then
        call this%xt3d%xt3d_ar(this%fmi%ibdgwfsat0, this%d11, this%id33, &
                               this%d33, this%fmi%gwfsat, this%id22, this%d22, &
                               this%iangle1, this%iangle2, this%iangle3, &
                               this%angle1, this%angle2, this%angle3)
      end if
    end if
    !
    ! -- Fill d11, d22, d33, angle1, angle2, angle3 using specific discharge
    call this%calcdispellipse()
    !
    ! -- Recalculate dispersion coefficients if the flows were updated
    if (this%fmi%iflowsupdated == 1) then
      if (this%ixt3d == 0) then
        call this%calcdispcoef()
      else if (this%ixt3d > 0) then
        call this%xt3d%xt3d_fcpc(this%dis%nodes, .false.)
      end if
    end if
  end subroutine dsp_ad

  !> @brief  Fill coefficient method for package
  !!
  !!  Method to calculate and fill coefficients for the package.
  !<
  subroutine dsp_fc(this, kiter, nodes, nja, matrix_sln, idxglo, rhs, cnew)
    ! -- modules
    ! -- dummy
    class(GwtDspType) :: this
    integer(I4B) :: kiter
    integer(I4B), intent(in) :: nodes
    integer(I4B), intent(in) :: nja
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(nja) :: idxglo
    real(DP), intent(inout), dimension(nodes) :: rhs
    real(DP), intent(inout), dimension(nodes) :: cnew
    ! -- local
    integer(I4B) :: n, m, idiag, idiagm, ipos, isympos, isymcon
    real(DP) :: dnm
    !
    if (this%ixt3d > 0) then
      call this%xt3d%xt3d_fc(kiter, matrix_sln, idxglo, rhs, cnew)
    else
      do n = 1, nodes
        if (this%fmi%ibdgwfsat0(n) == 0) cycle
        idiag = this%dis%con%ia(n)
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          if (this%dis%con%mask(ipos) == 0) cycle
          m = this%dis%con%ja(ipos)
          if (m < n) cycle
          if (this%fmi%ibdgwfsat0(m) == 0) cycle
          isympos = this%dis%con%jas(ipos)
          dnm = this%dispcoef(isympos)
          !
          ! -- Contribution to row n
          call matrix_sln%add_value_pos(idxglo(ipos), dnm)
          call matrix_sln%add_value_pos(idxglo(idiag), -dnm)
          !
          ! -- Contribution to row m
          idiagm = this%dis%con%ia(m)
          isymcon = this%dis%con%isym(ipos)
          call matrix_sln%add_value_pos(idxglo(isymcon), dnm)
          call matrix_sln%add_value_pos(idxglo(idiagm), -dnm)
        end do
      end do
    end if
  end subroutine dsp_fc

  !> @ brief Calculate flows for package
  !!
  !!  Method to calculate dispersion contribution to flowja
  !<
  subroutine dsp_cq(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(GwtDspType) :: this
    real(DP), intent(inout), dimension(:) :: cnew
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: n, m, ipos, isympos
    real(DP) :: dnm
    !
    ! -- Calculate dispersion and add to flowja
    if (this%ixt3d > 0) then
      call this%xt3d%xt3d_flowja(cnew, flowja)
    else
      do n = 1, this%dis%nodes
        if (this%fmi%ibdgwfsat0(n) == 0) cycle
        do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
          m = this%dis%con%ja(ipos)
          if (this%fmi%ibdgwfsat0(m) == 0) cycle
          isympos = this%dis%con%jas(ipos)
          dnm = this%dispcoef(isympos)
          flowja(ipos) = flowja(ipos) + dnm * (cnew(m) - cnew(n))
        end do
      end do
    end if
  end subroutine dsp_cq

  !> @ brief Allocate scalar variables for package
  !!
  !!  Method to allocate scalar variables for the package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GwtDspType) :: this
    ! -- local
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%idiffc, 'IDIFFC', this%memoryPath)
    call mem_allocate(this%idisp, 'IDISP', this%memoryPath)
    call mem_allocate(this%ialh, 'IALH', this%memoryPath)
    call mem_allocate(this%ialv, 'IALV', this%memoryPath)
    call mem_allocate(this%iath1, 'IATH1', this%memoryPath)
    call mem_allocate(this%iath2, 'IATH2', this%memoryPath)
    call mem_allocate(this%iatv, 'IATV', this%memoryPath)
    call mem_allocate(this%ixt3doff, 'IXT3DOFF', this%memoryPath)
    call mem_allocate(this%ixt3drhs, 'IXT3DRHS', this%memoryPath)
    call mem_allocate(this%ixt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%id22, 'ID22', this%memoryPath)
    call mem_allocate(this%id33, 'ID33', this%memoryPath)
    call mem_allocate(this%iangle1, 'IANGLE1', this%memoryPath)
    call mem_allocate(this%iangle2, 'IANGLE2', this%memoryPath)
    call mem_allocate(this%iangle3, 'IANGLE3', this%memoryPath)
    !
    ! -- Initialize
    this%idiffc = 0
    this%idisp = 0
    this%ialh = 0
    this%ialv = 0
    this%iath1 = 0
    this%iath2 = 0
    this%iatv = 0
    this%ixt3doff = 0
    this%ixt3drhs = 0
    this%ixt3d = 0
    this%id22 = 1
    this%id33 = 1
    this%iangle1 = 1
    this%iangle2 = 1
    this%iangle3 = 1
  end subroutine allocate_scalars

  !> @ brief Allocate arrays for package
  !!
  !!  Method to allocate arrays for the package.
  !<
  subroutine allocate_arrays(this, nodes)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GwtDspType) :: this
    integer(I4B), intent(in) :: nodes
    ! -- local
    !
    ! -- Allocate
    call mem_allocate(this%alh, nodes, 'ALH', trim(this%memoryPath))
    call mem_allocate(this%alv, nodes, 'ALV', trim(this%memoryPath))
    call mem_allocate(this%ath1, nodes, 'ATH1', trim(this%memoryPath))
    call mem_allocate(this%ath2, nodes, 'ATH2', trim(this%memoryPath))
    call mem_allocate(this%atv, nodes, 'ATV', trim(this%memoryPath))
    call mem_allocate(this%diffc, nodes, 'DIFFC', trim(this%memoryPath))
    call mem_allocate(this%d11, nodes, 'D11', trim(this%memoryPath))
    call mem_allocate(this%d22, nodes, 'D22', trim(this%memoryPath))
    call mem_allocate(this%d33, nodes, 'D33', trim(this%memoryPath))
    call mem_allocate(this%angle1, nodes, 'ANGLE1', trim(this%memoryPath))
    call mem_allocate(this%angle2, nodes, 'ANGLE2', trim(this%memoryPath))
    call mem_allocate(this%angle3, nodes, 'ANGLE3', trim(this%memoryPath))
    !
    ! -- Allocate dispersion coefficient array if xt3d not in use
    if (this%ixt3d == 0) then
      call mem_allocate(this%dispcoef, this%dis%njas, 'DISPCOEF', &
                        trim(this%memoryPath))
    else
      call mem_allocate(this%dispcoef, 0, 'DISPCOEF', trim(this%memoryPath))
    end if
  end subroutine allocate_arrays

  !> @ brief Deallocate memory
  !!
  !!  Method to deallocate memory for the package.
  !<
  subroutine dsp_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(GwtDspType) :: this
    ! -- local
    !
    ! -- Deallocate input memory
    call memorystore_remove(this%name_model, 'DSP', idm_context)
    !
    ! -- deallocate arrays
    if (this%inunit /= 0) then
      call mem_deallocate(this%alh, 'ALH', trim(this%memoryPath))
      call mem_deallocate(this%alv, 'ALV', trim(this%memoryPath))
      call mem_deallocate(this%ath1, 'ATH1', trim(this%memoryPath))
      call mem_deallocate(this%ath2, 'ATH2', trim(this%memoryPath))
      call mem_deallocate(this%atv, 'ATV', trim(this%memoryPath))
      call mem_deallocate(this%diffc)
      call mem_deallocate(this%d11)
      call mem_deallocate(this%d22)
      call mem_deallocate(this%d33)
      call mem_deallocate(this%angle1)
      call mem_deallocate(this%angle2)
      call mem_deallocate(this%angle3)
      call mem_deallocate(this%dispcoef)
      if (this%ixt3d > 0) call this%xt3d%xt3d_da()
    end if
    !
    ! -- deallocate objects
    if (this%ixt3d > 0) deallocate (this%xt3d)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%idiffc)
    call mem_deallocate(this%idisp)
    call mem_deallocate(this%ialh)
    call mem_deallocate(this%ialv)
    call mem_deallocate(this%iath1)
    call mem_deallocate(this%iath2)
    call mem_deallocate(this%iatv)
    call mem_deallocate(this%ixt3doff)
    call mem_deallocate(this%ixt3drhs)
    call mem_deallocate(this%ixt3d)
    call mem_deallocate(this%id22)
    call mem_deallocate(this%id33)
    call mem_deallocate(this%iangle1)
    call mem_deallocate(this%iangle2)
    call mem_deallocate(this%iangle3)
    !
    ! -- deallocate variables in NumericalPackageType
    call this%NumericalPackageType%da()
    !
    ! -- nullify pointers
    this%thetam => null()
  end subroutine dsp_da

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    use GwtDspInputModule, only: GwtDspParamFoundType
    class(GwTDspType) :: this
    type(GwtDspParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting DSP Options'
    write (this%iout, '(4x,a,i0)') 'XT3D formulation [0=INACTIVE, 1=ACTIVE, &
                                   &3=ACTIVE RHS] set to: ', this%ixt3d
    write (this%iout, '(1x,a,/)') 'End Setting DSP Options'
  end subroutine log_options

  !> @brief Update simulation mempath options
  !<
  subroutine source_options(this)
    ! -- modules
    !use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use GwtDspInputModule, only: GwtDspParamFoundType
    ! -- dummy
    class(GwtDspType) :: this
    ! -- locals
    type(GwtDspParamFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%ixt3doff, 'XT3D_OFF', this%input_mempath, &
                       found%xt3d_off)
    call mem_set_value(this%ixt3drhs, 'XT3D_RHS', this%input_mempath, &
                       found%xt3d_rhs)
    !
    ! -- set xt3d state flag
    if (found%xt3d_off) this%ixt3d = 0
    if (found%xt3d_rhs) this%ixt3d = 2
    !
    ! -- log options
    if (this%iout > 0) then
      call this%log_options(found)
    end if
  end subroutine source_options

  !> @brief Write dimensions to list file
  !<
  subroutine log_griddata(this, found)
    use GwtDspInputModule, only: GwtDspParamFoundType
    class(GwtDspType) :: this
    type(GwtDspParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting DSP Griddata'

    if (found%diffc) then
      write (this%iout, '(4x,a)') 'DIFFC set from input file'
    end if

    if (found%alh) then
      write (this%iout, '(4x,a)') 'ALH set from input file'
    end if

    if (found%alv) then
      write (this%iout, '(4x,a)') 'ALV set from input file'
    end if

    if (found%ath1) then
      write (this%iout, '(4x,a)') 'ATH1 set from input file'
    end if

    if (found%ath2) then
      write (this%iout, '(4x,a)') 'ATH2 set from input file'
    end if

    if (found%atv) then
      write (this%iout, '(4x,a)') 'ATV set from input file'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting DSP Griddata'

  end subroutine log_griddata

  !> @brief Update DSP simulation data from input mempath
  !<
  subroutine source_griddata(this)
    ! -- modules
    use SimModule, only: count_errors, store_error
    use MemoryManagerModule, only: mem_reallocate, mem_reassignptr
    use MemoryManagerExtModule, only: mem_set_value
    use ConstantsModule, only: LINELENGTH
    use GwtDspInputModule, only: GwtDspParamFoundType
    ! -- dummy
    class(GwtDsptype) :: this
    ! -- locals
    character(len=LINELENGTH) :: errmsg
    type(GwtDspParamFoundType) :: found
    integer(I4B), dimension(:), pointer, contiguous :: map
    ! -- formats
    !
    ! -- set map
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%diffc, 'DIFFC', this%input_mempath, map, found%diffc)
    call mem_set_value(this%alh, 'ALH', this%input_mempath, map, found%alh)
    call mem_set_value(this%alv, 'ALV', this%input_mempath, map, found%alv)
    call mem_set_value(this%ath1, 'ATH1', this%input_mempath, map, found%ath1)
    call mem_set_value(this%ath2, 'ATH2', this%input_mempath, map, found%ath2)
    call mem_set_value(this%atv, 'ATV', this%input_mempath, map, found%atv)
    !
    ! -- set active flags
    if (found%diffc) this%idiffc = 1
    if (found%alh) this%ialh = 1
    if (found%alv) this%ialv = 1
    if (found%ath1) this%iath1 = 1
    if (found%ath2) this%iath2 = 1
    if (found%atv) this%iatv = 1
    !
    ! -- reallocate diffc if not found
    if (.not. found%diffc) then
      call mem_reallocate(this%diffc, 0, 'DIFFC', trim(this%memoryPath))
    end if
    !
    ! -- set this%idisp flag
    if (found%alh) this%idisp = this%idisp + 1
    if (found%alv) this%idisp = this%idisp + 1
    if (found%ath1) this%idisp = this%idisp + 1
    if (found%ath2) this%idisp = this%idisp + 1
    !
    ! -- manage dispersion arrays
    if (this%idisp > 0) then
      if (.not. (found%alh .and. found%ath1)) then
        write (errmsg, '(a)') &
          'if dispersivities are specified then ALH and ATH1 are required.'
        call store_error(errmsg)
      end if
      ! -- If alv not specified then point it to alh
      if (.not. found%alv) &
        call mem_reassignptr(this%alv, 'ALV', trim(this%memoryPath), &
                             'ALH', trim(this%memoryPath))
      ! -- If ath2 not specified then point it to ath1
      if (.not. found%ath2) &
        call mem_reassignptr(this%ath2, 'ATH2', trim(this%memoryPath), &
                             'ATH1', trim(this%memoryPath))
      ! -- If atv not specified then point it to ath2
      if (.not. found%atv) &
        call mem_reassignptr(this%atv, 'ATV', trim(this%memoryPath), &
                             'ATH2', trim(this%memoryPath))
    else
      call mem_reallocate(this%alh, 0, 'ALH', trim(this%memoryPath))
      call mem_reallocate(this%alv, 0, 'ALV', trim(this%memoryPath))
      call mem_reallocate(this%ath1, 0, 'ATH1', trim(this%memoryPath))
      call mem_reallocate(this%ath2, 0, 'ATH2', trim(this%memoryPath))
      call mem_reallocate(this%atv, 0, 'ATV', trim(this%memoryPath))
    end if
    !
    ! -- log griddata
    if (this%iout > 0) then
      call this%log_griddata(found)
    end if
  end subroutine source_griddata

  !> @brief Calculate dispersion coefficients
  !<
  subroutine calcdispellipse(this)
    ! -- modules
    ! -- dummy
    class(GwtDspType) :: this
    ! -- local
    integer(I4B) :: nodes, n
    real(DP) :: q, qx, qy, qz
    real(DP) :: alh, alv, ath1, ath2, atv, a
    real(DP) :: al, at1, at2
    real(DP) :: qzoqsquared
    real(DP) :: dstar
    !
    ! -- loop through and calculate dispersion coefficients and angles
    nodes = size(this%d11)
    do n = 1, nodes
      !
      ! -- initialize
      this%d11(n) = DZERO
      this%d22(n) = DZERO
      this%d33(n) = DZERO
      this%angle1(n) = DZERO
      this%angle2(n) = DZERO
      this%angle3(n) = DZERO
      if (this%fmi%ibdgwfsat0(n) == 0) cycle
      !
      ! -- specific discharge
      qx = DZERO
      qy = DZERO
      qz = DZERO
      q = DZERO
      qx = this%fmi%gwfspdis(1, n)
      qy = this%fmi%gwfspdis(2, n)
      qz = this%fmi%gwfspdis(3, n)
      q = qx**2 + qy**2 + qz**2
      if (q > DZERO) q = sqrt(q)
      !
      ! -- dispersion coefficients
      alh = DZERO
      alv = DZERO
      ath1 = DZERO
      ath2 = DZERO
      atv = DZERO
      if (this%idisp > 0) then
        alh = this%alh(n)
        alv = this%alv(n)
        ath1 = this%ath1(n)
        ath2 = this%ath2(n)
        atv = this%atv(n)
      end if
      dstar = DZERO
      if (this%idiffc > 0) then
        ! -- Multiply diffusion coefficient by mobile porosity, defined
        !    as volume of voids in the mobile domain per aquifer volume
        dstar = this%diffc(n) * this%thetam(n)
      end if
      !
      ! -- Calculate the longitudal and transverse dispersivities
      al = DZERO
      at1 = DZERO
      at2 = DZERO
      if (q > DZERO) then
        qzoqsquared = (qz / q)**2
        al = alh * (DONE - qzoqsquared) + alv * qzoqsquared
        at1 = ath1 * (DONE - qzoqsquared) + atv * qzoqsquared
        at2 = ath2 * (DONE - qzoqsquared) + atv * qzoqsquared
      end if
      !
      ! -- Calculate and save the diagonal components of the dispersion tensor
      this%d11(n) = al * q + dstar
      this%d22(n) = at1 * q + dstar
      this%d33(n) = at2 * q + dstar
      !
      ! -- Angles of rotation if velocity based dispersion tensor
      if (this%idisp > 0) then
        !
        ! -- angles of rotation from model coordinates to direction of velocity
        ! qx / q = cos(a1) * cos(a2)
        ! qy / q = sin(a1) * cos(a2)
        ! qz / q = sin(a2)
        !
        ! -- angle3 is zero
        this%angle3(n) = DZERO
        !
        ! -- angle2
        a = DZERO
        if (q > DZERO) a = qz / q
        this%angle2(n) = asin(a)
        !
        ! -- angle1
        a = q * cos(this%angle2(n))
        if (a /= DZERO) then
          a = qx / a
        else
          a = DZERO
        end if
        !
        ! -- acos(1) not defined, so set to zero if necessary
        if (a <= -DONE) then
          this%angle1(n) = DPI
        elseif (a >= DONE) then
          this%angle1(n) = DZERO
        else
          this%angle1(n) = acos(a)
        end if
        !
      end if
    end do
  end subroutine calcdispellipse

  !> @brief Calculate dispersion coefficients
  !<
  subroutine calcdispcoef(this)
    ! -- modules
    use HGeoUtilModule, only: hyeff
    use GwfConductanceUtilsModule, only: staggered_thkfrac
    ! -- dummy
    class(GwtDspType) :: this
    ! -- local
    integer(I4B) :: nodes, n, m, idiag, ipos
    real(DP) :: clnm, clmn, dn, dm
    real(DP) :: vg1, vg2, vg3
    integer(I4B) :: ihc, isympos
    integer(I4B) :: iavgmeth
    real(DP) :: satn, satm, topn, topm, botn, botm
    real(DP) :: hwva, cond, cn, cm, denom
    real(DP) :: anm, amn, thksatn, thksatm
    !
    ! -- set iavgmeth = 1 to use arithmetic averaging for effective dispersion
    iavgmeth = 1
    !
    ! -- Process connections
    nodes = size(this%d11)
    do n = 1, nodes
      if (this%fmi%ibdgwfsat0(n) == 0) cycle
      idiag = this%dis%con%ia(n)
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        !
        ! -- Set m to connected cell
        m = this%dis%con%ja(ipos)
        !
        ! -- skip for lower triangle
        if (m < n) cycle
        isympos = this%dis%con%jas(ipos)
        this%dispcoef(isympos) = DZERO
        if (this%fmi%ibdgwfsat0(m) == 0) cycle
        !
        ! -- cell dimensions
        hwva = this%dis%con%hwva(isympos)
        clnm = this%dis%con%cl1(isympos)
        clmn = this%dis%con%cl2(isympos)
        ihc = this%dis%con%ihc(isympos)
        topn = this%dis%top(n)
        topm = this%dis%top(m)
        botn = this%dis%bot(n)
        botm = this%dis%bot(m)
        !
        ! -- flow model information
        satn = this%fmi%gwfsat(n)
        satm = this%fmi%gwfsat(m)
        !
        ! -- Calculate dispersion coefficient for cell n in the direction
        !    normal to the shared n-m face and for cell m in the direction
        !    normal to the shared n-m face.
        call this%dis%connection_normal(n, m, ihc, vg1, vg2, vg3, ipos)
        dn = hyeff(this%d11(n), this%d22(n), this%d33(n), &
                   this%angle1(n), this%angle2(n), this%angle3(n), &
                   vg1, vg2, vg3, iavgmeth)
        dm = hyeff(this%d11(m), this%d22(m), this%d33(m), &
                   this%angle1(m), this%angle2(m), this%angle3(m), &
                   vg1, vg2, vg3, iavgmeth)
        !
        ! -- Calculate dispersion conductance based on NPF subroutines and the
        !    effective dispersion coefficients dn and dm.
        if (ihc == 0) then
          clnm = satn * (topn - botn) * DHALF
          clmn = satm * (topm - botm) * DHALF
          anm = hwva
          !
          ! -- n is convertible and unsaturated
          if (satn == DZERO) then
            anm = DZERO
          else if (n > m .and. satn < DONE) then
            anm = DZERO
          end if
          !
          ! -- m is convertible and unsaturated
          if (satm == DZERO) then
            anm = DZERO
          else if (m > n .and. satm < DONE) then
            anm = DZERO
          end if
          !
          ! -- amn is the same as anm for vertical flow
          amn = anm
          !
        else
          !
          ! -- horizontal conductance
          !
          ! -- handle vertically staggered case
          if (ihc == 2) then
            thksatn = staggered_thkfrac(topn, botn, satn, topm, botm)
            thksatm = staggered_thkfrac(topm, botm, satm, topn, botn)
          else
            thksatn = (topn - botn) * satn
            thksatm = (topm - botm) * satm
          end if
          !
          ! -- calculate the saturated area term
          anm = thksatn * hwva
          amn = thksatm * hwva
          !
          ! -- n or m is unsaturated, so no dispersion
          if (satn == DZERO .or. satm == DZERO) then
            anm = DZERO
            amn = DZERO
          end if
          !
        end if
        !
        ! -- calculate conductance using the two half cell conductances
        cn = DZERO
        if (clnm > DZERO) cn = dn * anm / clnm
        cm = DZERO
        if (clmn > DZERO) cm = dm * amn / clmn
        denom = cn + cm
        if (denom > DZERO) then
          cond = cn * cm / denom
        else
          cond = DZERO
        end if
        !
        ! -- Assign the calculated dispersion conductance
        this%dispcoef(isympos) = cond
        !
      end do
    end do
  end subroutine calcdispcoef

end module GwtDspModule
