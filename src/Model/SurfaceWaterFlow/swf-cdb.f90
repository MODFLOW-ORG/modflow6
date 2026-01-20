!> @brief This module contains the CDB package methods
!!
!! This module can be used to represent outflow from streams using
!! a critical depth boundary.
!!
!<
module SwfCdbModule
  ! -- modules
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LENFTYPE, DNODATA, &
                             DPREC, DHALF, DTWO, DGRAVITY
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
  use SwfCxsModule, only: SwfCxsType
  !
  implicit none
  !
  private
  public :: cdb_create
  !
  character(len=LENFTYPE) :: ftype = 'CDB' !< package ftype
  character(len=16) :: text = '             CDB' !< package flow text string
  !
  type, extends(BndExtType) :: SwfCdbType

    integer(I4B), dimension(:), pointer, contiguous :: idcxs => null() !< cross section id
    real(DP), dimension(:), pointer, contiguous :: width => null() !< channel width
    real(DP), pointer :: gravconv => null() !< conversion factor gravity in m/s^2 to model units

    ! pointers other objects
    type(SwfCxsType), pointer :: cxs

  contains
    procedure :: allocate_scalars => cdb_allocate_scalars
    procedure :: allocate_arrays => cdb_allocate_arrays
    procedure :: source_options => cdb_options
    procedure :: log_cdb_options
    procedure :: bnd_rp => cdb_rp
    procedure :: bnd_cf => cdb_cf
    procedure :: bnd_fc => cdb_fc
    procedure :: bnd_da => cdb_da
    procedure :: define_listlabel
    procedure :: bound_value => cdb_bound_value
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => cdb_obs_supported
    procedure, public :: bnd_df_obs => cdb_df_obs
    procedure, public :: bnd_bd_obs => cdb_bd_obs
    ! -- methods for time series
    procedure, public :: bnd_rp_ts => cdb_rp_ts
    ! -- private
    procedure, private :: qcalc
  end type SwfCdbType

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new CDB Package object
  !!
  !<
  subroutine cdb_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath, dis, cxs, lengthconv, timeconv)
    ! -- dummy variables
    class(BndType), pointer :: packobj !< pointer to default package type
    integer(I4B), intent(in) :: id !< package id
    integer(I4B), intent(in) :: ibcnum !< boundary condition number
    integer(I4B), intent(in) :: inunit !< unit number of CDB package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< model name
    character(len=*), intent(in) :: pakname !< package name
    character(len=*), intent(in) :: mempath !< input mempath
    class(DisBaseType), pointer, intent(inout) :: dis !< the pointer to the discretization
    type(SwfCxsType), pointer, intent(in) :: cxs !< the pointer to the cxs package
    real(DP), intent(in) :: lengthconv !< conversion factor from model length to meters
    real(DP), intent(in) :: timeconv !< conversion factor from model time units to seconds
    ! -- local variables
    type(SwfCdbType), pointer :: cdbobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (cdbobj)
    packobj => cdbobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call cdbobj%allocate_scalars()
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

    ! -- store pointer to dis
    cdbobj%dis => dis

    ! -- store pointer to cxs
    cdbobj%cxs => cxs
    !
    ! -- store unit conversion
    cdbobj%gravconv = DGRAVITY * lengthconv * timeconv**2
  end subroutine cdb_create

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the CDB package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  subroutine cdb_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfcdbType object
    !
    ! -- call base type allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%gravconv, 'GRAVCONV', this%memoryPath)
    !
    ! -- Set values
    this%gravconv = DZERO
  end subroutine cdb_allocate_scalars

  !> @ brief Allocate arrays
    !!
    !! Allocate and initialize arrays for the SWF package
    !!
  !<
  subroutine cdb_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(SwfCdbType) :: this
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
    !
    ! -- checkin array input context pointer
    call mem_checkin(this%idcxs, 'IDCXS', this%memoryPath, &
                     'IDCXS', this%input_mempath)
    call mem_checkin(this%width, 'WIDTH', this%memoryPath, &
                     'WIDTH', this%input_mempath)
  end subroutine cdb_allocate_arrays

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate SWF package scalars and arrays.
  !!
  !<
  subroutine cdb_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfcdbType object
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- arrays
    call mem_deallocate(this%idcxs, 'IDCXS', this%memoryPath)
    call mem_deallocate(this%width, 'WIDTH', this%memoryPath)
    !
    ! -- scalars
    call mem_deallocate(this%gravconv)
  end subroutine cdb_da

  !> @ brief Source additional options for package
  !!
  !!  Source additional options for SWF package.
  !!
  !<
  subroutine cdb_options(this)
    ! -- modules
    use InputOutputModule, only: urword
    use MemoryManagerExtModule, only: mem_set_value
    use SwfCdbInputModule, only: SwfCdbParamFoundType
    ! -- dummy variables
    class(SwfCdbType), intent(inout) :: this !< SwfCdbType object
    ! -- local variables
    type(SwfCdbParamFoundType) :: found
    ! -- formats
    !
    ! -- source base BndExtType options
    call this%BndExtType%source_options()
    !
    ! -- source options from input context
    ! none
    !
    ! -- log SWF specific options
    call this%log_cdb_options(found)
  end subroutine cdb_options

  !> @ brief Log SWF specific package options
  !<
  subroutine log_cdb_options(this, found)
    ! -- modules
    use SwfCdbInputModule, only: SwfCdbParamFoundType
    ! -- dummy variables
    class(SwfCdbType), intent(inout) :: this !< BndExtType object
    type(SwfCdbParamFoundType), intent(in) :: found
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
  end subroutine log_cdb_options

  !> @ brief SWF read and prepare
    !!
  !<
  subroutine cdb_rp(this)
    use TdisModule, only: kper
    ! -- dummy
    class(SwfCdbType), intent(inout) :: this
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
  end subroutine cdb_rp

  !> @ brief Formulate the package hcof and rhs terms.
  !!
  !!  Formulate the hcof and rhs terms for the CDB package that will be
  !!  added to the coefficient matrix and right-hand side vector.
  !!
  !<
  subroutine cdb_cf(this)
    ! modules
    use MathUtilModule, only: get_perturbation
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfCdbType  object
    ! -- local variables
    integer(I4B) :: i, node
    real(DP) :: q
    real(DP) :: qeps
    real(DP) :: depth
    real(DP) :: derv
    real(DP) :: eps
    !
    ! -- Return if no inflows
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each cdb entry
    do i = 1, this%nbound

      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      ! -- calculate terms and add to hcof and rhs
      depth = this%xnew(node) - this%dis%bot(node)

      ! -- calculate unperturbed q
      q = this%qcalc(i, depth)

      ! -- calculate perturbed q
      eps = get_perturbation(depth)
      qeps = this%qcalc(i, depth + eps)

      ! -- calculate derivative
      derv = (qeps - q) / eps

      ! -- add terms to hcof and rhs
      this%hcof(i) = derv
      this%rhs(i) = -q + derv * this%xnew(node)

    end do
  end subroutine cdb_cf

  !> @brief Calculate critical depth boundary flow
  !<
  function qcalc(this, i, depth) result(q)
    ! modules
    ! dummy
    class(SwfCdbType) :: this
    integer(I4B), intent(in) :: i !< boundary number
    real(DP), intent(in) :: depth !< simulated depth (stage - elevation) in reach n for this iteration
    ! return
    real(DP) :: q
    ! local
    integer(I4B) :: idcxs
    real(DP) :: width
    real(DP) :: a
    real(DP) :: r

    idcxs = this%idcxs(i)
    width = this%width(i)
    a = this%cxs%get_area(idcxs, width, depth)
    r = this%cxs%get_hydraulic_radius(idcxs, width, depth, area=a)

    q = this%gravconv * a**DTWO * r
    if (q > DPREC) then
      q = q**DHALF
    else
      q = DZERO
    end if
    q = -q

  end function qcalc

  !> @ brief Copy hcof and rhs terms into solution.
  !!
  !!  Add the hcof and rhs terms for the CDB package to the
  !!  coefficient matrix and right-hand side vector.
  !!
  !<
  subroutine cdb_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfCdbType  object
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
      ! -- If mover is active and this cdb item is discharging,
      !    store available water (as positive value).
      if (this%imover == 1 .and. this%rhs(i) > DZERO) then
        call this%pakmvrobj%accumulate_qformvr(i, this%rhs(i))
      end if
    end do
  end subroutine cdb_fc

  !> @ brief Define the list label for the package
  !!
  !!  Method defined the list label for the CDB package. The list label is
  !!  the heading that is written to iout when PRINT_INPUT option is used.
  !!
  !<
  subroutine define_listlabel(this)
    ! -- dummy variables
    class(SwfCdbType), intent(inout) :: this !< SwfCdbType  object
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
  !! Function to determine if observations are supported by the CDB package.
  !! Observations are supported by the CDB package.
  !!
  !! @return  cdb_obs_supported       boolean indicating if observations are supported
  !!
  !<
  logical function cdb_obs_supported(this)
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfCdbType  object
    !
    ! -- set boolean
    cdb_obs_supported = .true.
  end function cdb_obs_supported

  !> @brief Define the observation types available in the package
  !!
  !! Method to define the observation types available in the CDB package.
  !!
  !<
  subroutine cdb_df_obs(this)
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfCdbType  object
    ! -- local variables
    integer(I4B) :: indx
    !
    ! -- initialize observations
    call this%obs%StoreObsType('cdb', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine cdb_df_obs

  !> @brief Save observations for the package
  !!
  !! Method to save simulated values for the CDB package.
  !!
  !<
  subroutine cdb_bd_obs(this)
    ! -- dummy variables
    class(SwfCdbType) :: this !< SwfCdbType  object
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
          case ('CDB')
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
  end subroutine cdb_bd_obs

  ! -- Procedure related to time series

  !> @brief Assign time series links for the package
  !!
  !! Assign the time series links for the CDB package. Only
  !! the Q variable can be defined with time series.
  !!
  !<
  subroutine cdb_rp_ts(this)
    ! -- dummy variables
    class(SwfCdbType), intent(inout) :: this !< SwfCdbType  object
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
  end subroutine cdb_rp_ts

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index
  !!  and row.
  !!
  !<
  function cdb_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(SwfCdbType), intent(inout) :: this !< BndExtType object
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
    case default
      errmsg = 'Programming error. CDB bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function cdb_bound_value

end module SwfCdbModule
