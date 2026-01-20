!> @brief This module contains the ZDG package methods
!!
!! This module can be used to represent outflow from streams using
!! a zero-depth-gradient boundary.
!!
!<
module SwfZdgModule
  ! -- modules
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LENFTYPE, DNODATA, DHALF
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use ObserveModule, only: ObserveType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType, &
                                  GetTimeSeriesLinkFromList
  use MatrixBaseModule
  use BaseDisModule, only: DisBaseType
  use Disv1dModule, only: Disv1dType
  use SwfCxsModule, only: SwfCxsType
  !
  implicit none
  !
  private
  public :: zdg_create
  !
  character(len=LENFTYPE) :: ftype = 'ZDG' !< package ftype
  character(len=16) :: text = '             ZDG' !< package flow text string
  !
  type, extends(BndExtType) :: SwfZdgType

    integer(I4B), dimension(:), pointer, contiguous :: idcxs => null() !< cross section id
    real(DP), dimension(:), pointer, contiguous :: width => null() !< channel width
    real(DP), dimension(:), pointer, contiguous :: slope => null() !< channel slope
    real(DP), dimension(:), pointer, contiguous :: rough => null() !< channel roughness
    real(DP), pointer :: unitconv => null() !< conversion factor for roughness to length and time units of meters and seconds

    ! -- pointers other objects
    type(Disv1dType), pointer :: disv1d
    type(SwfCxsType), pointer :: cxs

  contains
    procedure :: allocate_scalars => zdg_allocate_scalars
    procedure :: allocate_arrays => zdg_allocate_arrays
    procedure :: source_options => zdg_options
    procedure :: log_zdg_options
    procedure :: bnd_rp => zdg_rp
    procedure :: bnd_cf => zdg_cf
    procedure :: bnd_fc => zdg_fc
    procedure :: bnd_da => zdg_da
    procedure :: define_listlabel
    procedure :: bound_value => zdg_bound_value
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => zdg_obs_supported
    procedure, public :: bnd_df_obs => zdg_df_obs
    procedure, public :: bnd_bd_obs => zdg_bd_obs
    ! -- methods for time series
    procedure, public :: bnd_rp_ts => zdg_rp_ts
    ! -- private
    procedure, private :: qcalc
  end type SwfZdgType

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new ZDG Package object
  !!
  !<
  subroutine zdg_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath, dis, cxs, unitconv)
    ! -- dummy variables
    class(BndType), pointer :: packobj !< pointer to default package type
    integer(I4B), intent(in) :: id !< package id
    integer(I4B), intent(in) :: ibcnum !< boundary condition number
    integer(I4B), intent(in) :: inunit !< unit number of ZDG package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< model name
    character(len=*), intent(in) :: pakname !< package name
    character(len=*), intent(in) :: mempath !< input mempath
    class(DisBaseType), pointer, intent(inout) :: dis !< the pointer to the discretization
    type(SwfCxsType), pointer, intent(in) :: cxs !< the pointer to the cxs package
    real(DP), intent(in) :: unitconv !< unit conversion for roughness
    ! -- local variables
    type(SwfZdgType), pointer :: zdgobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (zdgobj)
    packobj => zdgobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call zdgobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1
    packobj%ictMemPath = create_mem_path(namemodel, 'DFW')
    !
    ! -- store pointer to disv1d
    select type (dis)
    type is (Disv1dType)
      zdgobj%disv1d => dis
    end select
    !
    ! -- store pointer to cxs
    zdgobj%cxs => cxs
    !
    ! -- store unit conversion
    zdgobj%unitconv = unitconv
  end subroutine zdg_create

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the ZDG package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  subroutine zdg_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(SwfZdgType) :: this !< SwfZdgType object
    !
    ! -- call base type allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%unitconv, 'UNITCONV', this%memoryPath)
    !
    ! -- Set values
    this%unitconv = DZERO
  end subroutine zdg_allocate_scalars

  !> @ brief Allocate arrays
    !!
    !! Allocate and initialize arrays for the SWF package
    !!
  !<
  subroutine zdg_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(SwfZdgType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    !
    ! -- call BndExtType allocate scalars
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set array input context pointer
    call mem_setptr(this%idcxs, 'IDCXS', this%input_mempath)
    call mem_setptr(this%width, 'WIDTH', this%input_mempath)
    call mem_setptr(this%slope, 'SLOPE', this%input_mempath)
    call mem_setptr(this%rough, 'ROUGH', this%input_mempath)
    !
    ! -- checkin array input context pointer
    call mem_checkin(this%idcxs, 'IDCXS', this%memoryPath, &
                     'IDCXS', this%input_mempath)
    call mem_checkin(this%width, 'WIDTH', this%memoryPath, &
                     'WIDTH', this%input_mempath)
    call mem_checkin(this%slope, 'SLOPE', this%memoryPath, &
                     'SLOPE', this%input_mempath)
    call mem_checkin(this%rough, 'ROUGH', this%memoryPath, &
                     'ROUGH', this%input_mempath)
  end subroutine zdg_allocate_arrays

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate SWF package scalars and arrays.
  !!
  !<
  subroutine zdg_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(SwfZdgType) :: this !< SwfZdgType object
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- arrays
    call mem_deallocate(this%idcxs, 'IDCXS', this%memoryPath)
    call mem_deallocate(this%width, 'WIDTH', this%memoryPath)
    call mem_deallocate(this%slope, 'SLOPE', this%memoryPath)
    call mem_deallocate(this%rough, 'ROUGH', this%memoryPath)
    !
    ! -- scalars
    call mem_deallocate(this%unitconv)
  end subroutine zdg_da

  !> @ brief Source additional options for package
  !!
  !!  Source additional options for SWF package.
  !!
  !<
  subroutine zdg_options(this)
    ! -- modules
    use InputOutputModule, only: urword
    use MemoryManagerExtModule, only: mem_set_value
    use SwfZdgInputModule, only: SwfZdgParamFoundType
    ! -- dummy variables
    class(SwfZdgType), intent(inout) :: this !< SwfZdgType object
    ! -- local variables
    type(SwfZdgParamFoundType) :: found
    ! -- formats
    !
    ! -- source base BndExtType options
    call this%BndExtType%source_options()
    !
    ! -- source options from input context
    ! none
    !
    ! -- log SWF specific options
    call this%log_zdg_options(found)
  end subroutine zdg_options

  !> @ brief Log SWF specific package options
  !<
  subroutine log_zdg_options(this, found)
    ! -- modules
    use SwfZdgInputModule, only: SwfZdgParamFoundType
    ! -- dummy variables
    class(SwfZdgType), intent(inout) :: this !< BndExtType object
    type(SwfZdgParamFoundType), intent(in) :: found
    ! -- local variables
    ! -- format
    !
    ! -- log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'
    !
    ! if (found%mover) then
    !   write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
    ! end if
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine log_zdg_options

  !> @ brief SWF read and prepare
    !!
  !<
  subroutine zdg_rp(this)
    use TdisModule, only: kper
    ! -- dummy
    class(SwfZdgType), intent(inout) :: this
    ! -- local
    !
    if (this%iper /= kper) return
    !
    ! -- Call the parent class read and prepare
    call this%BndExtType%bnd_rp()
    !
    ! -- Write the list to iout if requested
    if (this%iprpak /= 0) then
      call this%write_list()
    end if
  end subroutine zdg_rp

  !> @ brief Formulate the package hcof and rhs terms.
  !!
  !!  Formulate the hcof and rhs terms for the ZDG package that will be
  !!  added to the coefficient matrix and right-hand side vector.
  !!
  !<
  subroutine zdg_cf(this)
    ! modules
    use MathUtilModule, only: get_perturbation
    use SmoothingModule, only: sQuadratic
    ! dummy variables
    class(SwfZdgType) :: this !< SwfZdgType  object
    ! local variables
    integer(I4B) :: i, node
    real(DP) :: q
    real(DP) :: qeps
    real(DP) :: absdhdxsq
    real(DP) :: depth
    real(DP) :: derv
    real(DP) :: eps
    real(DP) :: range = 1.d-6
    real(DP) :: dydx
    real(DP) :: smooth_factor
    !
    ! -- Return if no inflows
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each zdg entry
    do i = 1, this%nbound

      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      ! -- calculate terms and add to hcof and rhs
      absdhdxsq = this%slope(i)**DHALF
      depth = this%xnew(node) - this%dis%bot(node)

      ! smooth the depth
      call sQuadratic(depth, range, dydx, smooth_factor)
      depth = depth * smooth_factor

      ! -- calculate unperturbed q
      q = -this%qcalc(i, depth, this%unitconv)

      ! -- calculate perturbed q
      eps = get_perturbation(depth)
      qeps = -this%qcalc(i, depth + eps, this%unitconv)

      ! -- calculate derivative
      derv = (qeps - q) / eps

      ! -- add terms to hcof and rhs
      this%hcof(i) = derv
      this%rhs(i) = -q + derv * this%xnew(node)

    end do
  end subroutine zdg_cf

  ! !> @brief Calculate flow
  !!
  !! Calculate volumetric flow rate for the zero-depth gradient
  !! condition.  Flow is positive.
  ! !<
  function qcalc(this, i, depth, unitconv) result(q)
    ! dummy
    class(SwfZdgType) :: this
    integer(I4B), intent(in) :: i !< boundary number
    real(DP), intent(in) :: depth !< simulated depth (stage - elevation) in reach n for this iteration
    real(DP), intent(in) :: unitconv !< conversion factor for roughness to length and time units of meters and seconds
    ! return
    real(DP) :: q
    ! local
    integer(I4B) :: idcxs
    real(DP) :: width
    real(DP) :: rough
    real(DP) :: slope
    real(DP) :: conveyance

    idcxs = this%idcxs(i)
    width = this%width(i)
    rough = this%rough(i)
    slope = this%slope(i)
    conveyance = this%cxs%get_conveyance(idcxs, width, depth, rough)
    q = conveyance * slope**DHALF * unitconv

  end function qcalc

  !> @ brief Copy hcof and rhs terms into solution.
  !!
  !!  Add the hcof and rhs terms for the ZDG package to the
  !!  coefficient matrix and right-hand side vector.
  !!
  !<
  subroutine zdg_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy variables
    class(SwfZdgType) :: this !< SwfZdgType  object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! -- local variables
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: ipos
    !
    ! -- pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
      !
      ! -- If mover is active and this zdg item is discharging,
      !    store available water (as positive value).
      if (this%imover == 1 .and. this%rhs(i) > DZERO) then
        call this%pakmvrobj%accumulate_qformvr(i, this%rhs(i))
      end if
    end do
  end subroutine zdg_fc

  !> @ brief Define the list label for the package
  !!
  !!  Method defined the list label for the ZDG package. The list label is
  !!  the heading that is written to iout when PRINT_INPUT option is used.
  !!
  !<
  subroutine define_listlabel(this)
    ! -- dummy variables
    class(SwfZdgType), intent(inout) :: this !< SwfZdgType  object
    !
    ! -- create the header list label
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'FLOW RATE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  ! -- Procedures related to observations

  !> @brief Determine if observations are supported.
  !!
  !! Function to determine if observations are supported by the ZDG package.
  !! Observations are supported by the ZDG package.
  !!
  !! @return  zdg_obs_supported       boolean indicating if observations are supported
  !!
  !<
  logical function zdg_obs_supported(this)
    ! -- dummy variables
    class(SwfZdgType) :: this !< SwfZdgType  object
    !
    ! -- set boolean
    zdg_obs_supported = .true.
  end function zdg_obs_supported

  !> @brief Define the observation types available in the package
  !!
  !! Method to define the observation types available in the ZDG package.
  !!
  !<
  subroutine zdg_df_obs(this)
    ! -- dummy variables
    class(SwfZdgType) :: this !< SwfZdgType  object
    ! -- local variables
    integer(I4B) :: indx
    !
    ! -- initialize observations
    call this%obs%StoreObsType('zdg', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine zdg_df_obs

  !> @brief Save observations for the package
  !!
  !! Method to save simulated values for the ZDG package.
  !!
  !<
  subroutine zdg_bd_obs(this)
    ! -- dummy variables
    class(SwfZdgType) :: this !< SwfZdgType  object
    ! -- local variables
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: jj
    real(DP) :: v
    type(ObserveType), pointer :: obsrv => null()
    !
    ! -- clear the observations
    call this%obs%obs_bd_clear()
    !
    ! -- Save simulated values for all of package's observations.
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      if (obsrv%BndFound) then
        do n = 1, obsrv%indxbnds_count
          v = DNODATA
          jj = obsrv%indxbnds(n)
          select case (obsrv%ObsTypeId)
          case ('TO-MVR')
            if (this%imover == 1) then
              v = this%pakmvrobj%get_qtomvr(jj)
              if (v > DZERO) then
                v = -v
              end if
            end if
          case ('ZDG')
            v = this%simvals(jj)
          case default
            errmsg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
            call store_error(errmsg)
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      else
        call this%obs%SaveOneSimval(obsrv, DNODATA)
      end if
    end do
  end subroutine zdg_bd_obs

  ! -- Procedure related to time series

  !> @brief Assign time series links for the package
  !!
  !! Assign the time series links for the ZDG package. Only
  !! the Q variable can be defined with time series.
  !!
  !<
  subroutine zdg_rp_ts(this)
    ! -- dummy variables
    class(SwfZdgType), intent(inout) :: this !< SwfZdgType  object
    ! -- local variables
    integer(I4B) :: i, nlinks
    type(TimeSeriesLinkType), pointer :: tslink => null()
    !
    ! -- set up the time series links
    nlinks = this%TsManager%boundtslinks%Count()
    do i = 1, nlinks
      tslink => GetTimeSeriesLinkFromList(this%TsManager%boundtslinks, i)
      if (associated(tslink)) then
        if (tslink%JCol == 1) then
          tslink%Text = 'Q'
        end if
      end if
    end do
  end subroutine zdg_rp_ts

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index
  !!  and row.
  !!
  !<
  function zdg_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(SwfZdgType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%idcxs(row)
    case (2)
      bndval = this%width(row)
    case (3)
      bndval = this%slope(row)
    case (4)
      bndval = this%rough(row)
    case default
      errmsg = 'Programming error. ZDG bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function zdg_bound_value

end module SwfZdgModule
