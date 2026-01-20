!> @brief This module contains the FLW package methods
!!
!! This module can be used to represent inflow to streams.  It is
!! designed similarly to the GWF WEL package.
!!
!<
module SwfFlwModule
  ! -- modules
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LENFTYPE, DNODATA
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use ObserveModule, only: ObserveType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType, &
                                  GetTimeSeriesLinkFromList
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: flw_create
  !
  character(len=LENFTYPE) :: ftype = 'FLW' !< package ftype
  character(len=16) :: text = '             FLW' !< package flow text string
  !
  type, extends(BndExtType) :: SwfFlwType
    real(DP), dimension(:), pointer, contiguous :: q => null() !< volumetric rate
  contains
    procedure :: allocate_scalars => flw_allocate_scalars
    procedure :: allocate_arrays => flw_allocate_arrays
    procedure :: source_options => flw_options
    procedure :: log_flw_options
    procedure :: bnd_rp => flw_rp
    procedure :: bnd_cf => flw_cf
    procedure :: bnd_fc => flw_fc
    procedure :: bnd_da => flw_da
    procedure :: define_listlabel
    procedure :: bound_value => flw_bound_value
    procedure :: q_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => flw_obs_supported
    procedure, public :: bnd_df_obs => flw_df_obs
    procedure, public :: bnd_bd_obs => flw_bd_obs
    ! -- methods for time series
    procedure, public :: bnd_rp_ts => flw_rp_ts
  end type SwfFlwType

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new FLW Package object
  !!
  !<
  subroutine flw_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath)
    ! -- dummy variables
    class(BndType), pointer :: packobj !< pointer to default package type
    integer(I4B), intent(in) :: id !< package id
    integer(I4B), intent(in) :: ibcnum !< boundary condition number
    integer(I4B), intent(in) :: inunit !< unit number of FLW package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< model name
    character(len=*), intent(in) :: pakname !< package name
    character(len=*), intent(in) :: mempath !< input mempath
    ! -- local variables
    type(SwfFlwType), pointer :: flwobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (flwobj)
    packobj => flwobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call flwobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1
    packobj%ictMemPath = ''
  end subroutine flw_create

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the FLW package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  subroutine flw_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType object
    !
    ! -- call base type allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    ! none for this package
    !
    ! -- Set values
    ! none for this package
  end subroutine flw_allocate_scalars

  !> @ brief Allocate arrays
    !!
    !! Allocate and initialize arrays for the SWF package
    !!
  !<
  subroutine flw_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(SwfFlwType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    !
    ! -- call BndExtType allocate scalars
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set q array input context pointer
    call mem_setptr(this%q, 'Q', this%input_mempath)
    !
    ! -- checkin q array input context pointer
    call mem_checkin(this%q, 'Q', this%memoryPath, &
                     'Q', this%input_mempath)
  end subroutine flw_allocate_arrays

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate SWF package scalars and arrays.
  !!
  !<
  subroutine flw_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType object
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- scalars
    call mem_deallocate(this%q, 'Q', this%memoryPath)
  end subroutine flw_da

  !> @ brief Source additional options for package
  !!
  !!  Source additional options for SWF package.
  !!
  !<
  subroutine flw_options(this)
    ! -- modules
    use InputOutputModule, only: urword
    use MemoryManagerExtModule, only: mem_set_value
    use SwfFlwInputModule, only: SwfFlwParamFoundType
    ! -- dummy variables
    class(SwfFlwType), intent(inout) :: this !< SwfFlwType object
    ! -- local variables
    type(SwfFlwParamFoundType) :: found
    ! -- formats
    !
    ! -- source base BndExtType options
    call this%BndExtType%source_options()
    !
    ! -- source options from input context
    ! none
    !
    ! -- log SWF specific options
    call this%log_flw_options(found)
  end subroutine flw_options

  !> @ brief Log SWF specific package options
  !<
  subroutine log_flw_options(this, found)
    ! -- modules
    use SwfFlwInputModule, only: SwfFlwParamFoundType
    ! -- dummy variables
    class(SwfFlwType), intent(inout) :: this !< BndExtType object
    type(SwfFlwParamFoundType), intent(in) :: found
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
  end subroutine log_flw_options

  !> @ brief SWF read and prepare
    !!
  !<
  subroutine flw_rp(this)
    use TdisModule, only: kper
    ! -- dummy
    class(SwfFlwType), intent(inout) :: this
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
  end subroutine flw_rp

  !> @ brief Formulate the package hcof and rhs terms.
  !!
  !!  Formulate the hcof and rhs terms for the FLW package that will be
  !!  added to the coefficient matrix and right-hand side vector.
  !!
  !<
  subroutine flw_cf(this)
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType  object
    ! -- local variables
    integer(I4B) :: i, node
    real(DP) :: q
    !
    ! -- Return if no inflows
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each flw entry
    do i = 1, this%nbound
      node = this%nodelist(i)
      this%hcof(i) = DZERO
      if (this%ibound(node) <= 0) then
        this%rhs(i) = DZERO
        cycle
      end if
      q = this%q_mult(i)
      this%rhs(i) = -q
    end do
  end subroutine flw_cf

  !> @ brief Copy hcof and rhs terms into solution.
  !!
  !!  Add the hcof and rhs terms for the FLW package to the
  !!  coefficient matrix and right-hand side vector.
  !!
  !<
  subroutine flw_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType  object
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
      ! -- If mover is active and this flw item is discharging,
      !    store available water (as positive value).
      if (this%imover == 1 .and. this%rhs(i) > DZERO) then
        call this%pakmvrobj%accumulate_qformvr(i, this%rhs(i))
      end if
    end do
  end subroutine flw_fc

  !> @ brief Define the list label for the package
  !!
  !!  Method defined the list label for the FLW package. The list label is
  !!  the heading that is written to iout when PRINT_INPUT option is used.
  !!
  !<
  subroutine define_listlabel(this)
    ! -- dummy variables
    class(SwfFlwType), intent(inout) :: this !< SwfFlwType  object
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
  !! Function to determine if observations are supported by the FLW package.
  !! Observations are supported by the FLW package.
  !!
  !! @return  flw_obs_supported       boolean indicating if observations are supported
  !!
  !<
  logical function flw_obs_supported(this)
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType  object
    !
    ! -- set boolean
    flw_obs_supported = .true.
  end function flw_obs_supported

  !> @brief Define the observation types available in the package
  !!
  !! Method to define the observation types available in the FLW package.
  !!
  !<
  subroutine flw_df_obs(this)
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType  object
    ! -- local variables
    integer(I4B) :: indx
    !
    ! -- initialize observations
    call this%obs%StoreObsType('flw', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine flw_df_obs

  !> @brief Save observations for the package
  !!
  !! Method to save simulated values for the FLW package.
  !!
  !<
  subroutine flw_bd_obs(this)
    ! -- dummy variables
    class(SwfFlwType) :: this !< SwfFlwType  object
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
          case ('FLW')
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
  end subroutine flw_bd_obs

  ! -- Procedure related to time series

  !> @brief Assign time series links for the package
  !!
  !! Assign the time series links for the FLW package. Only
  !! the Q variable can be defined with time series.
  !!
  !<
  subroutine flw_rp_ts(this)
    ! -- dummy variables
    class(SwfFlwType), intent(inout) :: this !< SwfFlwType  object
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
  end subroutine flw_rp_ts

  function q_mult(this, row) result(q)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(SwfFlwType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: q
    !
    if (this%iauxmultcol > 0) then
      q = this%q(row) * this%auxvar(this%iauxmultcol, row)
    else
      q = this%q(row)
    end if
  end function q_mult

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index
  !!  and row.
  !!
  !<
  function flw_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(SwfFlwType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%q_mult(row)
    case default
      errmsg = 'Programming error. FLW bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function flw_bound_value

end module SwfFlwModule
