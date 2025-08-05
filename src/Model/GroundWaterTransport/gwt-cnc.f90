module GwtCncModule
  !
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, NAMEDBOUNDFLAG, LENFTYPE, &
                             LENPACKAGENAME, LENVARNAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use ObsModule, only: DefaultObsIdProcessor
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObserveModule, only: ObserveType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType, &
                                  GetTimeSeriesLinkFromList
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: cnc_create
  !
  character(len=LENFTYPE) :: ftype = 'CNC'
  character(len=LENPACKAGENAME) :: text = '             CNC'
  !
  type, extends(BndExtType) :: GwtCncType

    real(DP), dimension(:), pointer, contiguous :: tspvar => null() !< constant concentration array
    real(DP), dimension(:), pointer, contiguous :: ratecncin => null() !simulated flows into constant conc (excluding other concs)
    real(DP), dimension(:), pointer, contiguous :: ratecncout => null() !simulated flows out of constant conc (excluding to other concs)
    character(len=LENVARNAME) :: depvartype = '' !< stores string of dependent variable type, depending on model type
  contains
    procedure :: bnd_rp => cnc_rp
    procedure :: bnd_ad => cnc_ad
    procedure :: bnd_ck => cnc_ck
    procedure :: bnd_fc => cnc_fc
    procedure :: bnd_cq => cnc_cq
    procedure :: bnd_bd => cnc_bd
    procedure :: bnd_da => cnc_da
    procedure :: allocate_arrays => cnc_allocate_arrays
    procedure :: define_listlabel
    procedure :: bound_value => cnc_bound_value
    procedure :: conc_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => cnc_obs_supported
    procedure, public :: bnd_df_obs => cnc_df_obs
    ! -- method for time series
    procedure, public :: bnd_rp_ts => cnc_rp_ts
  end type GwtCncType

contains

  !> @brief Create a new constant concentration or temperature package
  !!
  !! Routine points packobj to the newly created package
  !<
  subroutine cnc_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        depvartype, mempath)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=LENVARNAME), intent(in) :: depvartype
    character(len=*), intent(in) :: mempath
    ! -- local
    type(GwtCncType), pointer :: cncobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (cncobj)
    packobj => cncobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call cncobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    ! -- store values
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    !
    ! -- Store the appropriate label based on the dependent variable
    cncobj%depvartype = depvartype
  end subroutine cnc_create

  !> @brief Allocate arrays specific to the constant concentration/tempeature
  !! package.
  !<
  subroutine cnc_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(GwtCncType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    integer(I4B) :: i
    !
    ! -- call standard BndType allocate scalars
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- allocate ratecncex
    call mem_allocate(this%ratecncin, this%maxbound, 'RATECNCIN', this%memoryPath)
    call mem_allocate(this%ratecncout, this%maxbound, 'RATECNCOUT', &
                      this%memoryPath)
    do i = 1, this%maxbound
      this%ratecncin(i) = DZERO
      this%ratecncout(i) = DZERO
    end do

    ! -- set input context pointers
    call mem_setptr(this%tspvar, 'TSPVAR', this%input_mempath)

    ! -- checkin input context pointers
    call mem_checkin(this%tspvar, 'TSPVAR', this%memoryPath, &
                     'TSPVAR', this%input_mempath)
  end subroutine cnc_allocate_arrays

  !> @brief Constant concentration/temperature read and prepare (rp) routine
  !<
  subroutine cnc_rp(this)
    ! -- modules
    use SimModule, only: store_error
    use InputOutputModule, only: lowcase
    implicit none
    ! -- dummy
    class(GwtCncType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, node, ibd, ierr
    character(len=30) :: nodestr
    character(len=LENVARNAME) :: dvtype
    !
    ! -- Reset previous CNCs to active cell
    do i = 1, this%nbound
      node = this%nodelist(i)
      this%ibound(node) = this%ibcnum
    end do
    !
    ! -- Call the parent class read and prepare
    call this%BndExtType%bnd_rp()
    !
    ! -- Set ibound to -(ibcnum + 1) for constant concentration cells
    ierr = 0
    do i = 1, this%nbound
      node = this%nodelist(i)
      ibd = this%ibound(node)
      if (ibd < 0) then
        call this%dis%noder_to_string(node, nodestr)
        dvtype = trim(this%depvartype)
        call lowcase(dvtype)
        call store_error('Cell is already a constant ' &
                         //dvtype//': '//trim(adjustl(nodestr)))
        ierr = ierr + 1
      else
        this%ibound(node) = -this%ibcnum
      end if
    end do
    !
    ! -- Stop if errors detected
    if (ierr > 0) then
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Write the list to iout if requested
    if (this%iprpak /= 0) then
      call this%write_list()
    end if
  end subroutine cnc_rp

  !> @brief Constant concentration/temperature package advance routine
  !!
  !! Add package connections to matrix
  !<
  subroutine cnc_ad(this)
    ! -- modules
    ! -- dummy
    class(GwtCncType) :: this
    ! -- local
    integer(I4B) :: i, node
    real(DP) :: cb
    ! -- formats
    !
    ! -- Advance the time series
    call this%TsManager%ad()
    !
    ! -- Process each entry in the constant concentration cell list
    do i = 1, this%nbound
      node = this%nodelist(i)
      cb = this%conc_mult(i)
      !
      this%xnew(node) = cb
      this%xold(node) = this%xnew(node)
    end do
    !
    ! -- For each observation, push simulated value and corresponding
    !    simulation time from "current" to "preceding" and reset
    !    "current" value.
    call this%obs%obs_ad()
  end subroutine cnc_ad

  !> @brief Check constant concentration/temperature boundary condition data
  !<
  subroutine cnc_ck(this)
    ! -- modules
    ! -- dummy
    class(GwtCncType), intent(inout) :: this
    ! -- local
    character(len=30) :: nodestr
    integer(I4B) :: i
    integer(I4B) :: node
    ! -- formats
    character(len=*), parameter :: fmtcncerr = &
      &"('Specified dependent variable boundary ',i0, &
      &' conc (',g0,') is less than zero for cell', a)"
    !
    ! -- check stress period data
    do i = 1, this%nbound
      node = this%nodelist(i)
      ! -- accumulate errors
      if (this%conc_mult(i) < DZERO) then
        call this%dis%noder_to_string(node, nodestr)
        write (errmsg, fmt=fmtcncerr) i, this%tspvar(i), trim(nodestr)
        call store_error(errmsg)
      end if
    end do
    !
    ! -- write summary of cnc package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine cnc_ck

  !> @brief Override bnd_fc and do nothing
  !!
  !! For constant concentration/temperature boundary type, the call to bnd_fc
  !! needs to be overwritten to prevent logic found in bnd from being executed
  !<
  subroutine cnc_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GwtCncType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
  end subroutine cnc_fc

  !> @brief Calculate flow associated with constant concentration/temperature
  !! boundary
  !!
  !! This method overrides bnd_cq()
  !<
  subroutine cnc_cq(this, x, flowja, iadv)
    ! -- modules
    ! -- dummy
    class(GwtCncType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: ipos
    integer(I4B) :: node
    integer(I4B) :: n2
    integer(I4B) :: idiag
    real(DP) :: rate
    real(DP) :: ratein, rateout
    real(DP) :: q
    !
    ! -- If no boundaries, skip flow calculations.
    if (this%nbound > 0) then
      !
      ! -- Loop through each boundary calculating flow.
      do i = 1, this%nbound
        node = this%nodelist(i)
        idiag = this%dis%con%ia(node)
        rate = DZERO
        ratein = DZERO
        rateout = DZERO
        !
        ! -- Calculate the flow rate into the cell.
        do ipos = this%dis%con%ia(node) + 1, &
          this%dis%con%ia(node + 1) - 1
          q = flowja(ipos)
          rate = rate - q
          ! -- only accumulate chin and chout for active
          !    connected cells
          n2 = this%dis%con%ja(ipos)
          if (this%ibound(n2) > 0) then
            if (q < DZERO) then
              ratein = ratein - q
            else
              rateout = rateout + q
            end if
          end if
        end do
        !
        ! -- For CNC, store total flow in rhs so it is available for other
        !    calculations
        this%rhs(i) = -rate
        this%hcof(i) = DZERO
        !
        ! -- Save simulated value to simvals array.
        this%simvals(i) = rate
        this%ratecncin(i) = ratein
        this%ratecncout(i) = rateout
        flowja(idiag) = flowja(idiag) + rate
        !
      end do
      !
    end if
  end subroutine cnc_cq

  !> @brief Add package ratin/ratout to model budget
  !<
  subroutine cnc_bd(this, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(GwtCncType) :: this
    ! -- local
    type(BudgetType), intent(inout) :: model_budget
    real(DP) :: ratin
    real(DP) :: ratout
    real(DP) :: dum
    integer(I4B) :: isuppress_output
    !
    isuppress_output = 0
    call rate_accumulator(this%ratecncin(1:this%nbound), ratin, dum)
    call rate_accumulator(this%ratecncout(1:this%nbound), ratout, dum)
    call model_budget%addentry(ratin, ratout, delt, this%text, &
                               isuppress_output, this%packName)
  end subroutine cnc_bd

  !> @brief Deallocate memory
  !!
  !!  Method to deallocate memory for the package.
  !<
  subroutine cnc_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwtCncType) :: this
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- arrays
    call mem_deallocate(this%ratecncin)
    call mem_deallocate(this%ratecncout)
    call mem_deallocate(this%tspvar, 'TSPVAR', this%memoryPath)
  end subroutine cnc_da

  !> @brief Define labels used in list file
  !!
  !! Define the list heading that is written to iout when PRINT_INPUT option
  !! is used.
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(GwtCncType), intent(inout) :: this
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), &
      trim(this%depvartype)
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  !> @brief Procedure related to observation processing
  !!
  !! This routine:
  !!   - returns true because the SDV package supports observations,
  !!   - overrides packagetype%_obs_supported()
  logical function cnc_obs_supported(this)
    ! -- dummy
    class(GwtCncType) :: this
    !
    cnc_obs_supported = .true.
  end function cnc_obs_supported

  !> @brief Procedure related to observation processing
  !!
  !! This routine:
  !!   - defines observations
  !!   - stores observation types supported by either of the SDV packages
  !!     (CNC or CNT),
  !!   - overrides BndExtType%bnd_df_obs
  !<
  subroutine cnc_df_obs(this)
    ! -- dummy
    class(GwtCncType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType(this%filtyp, .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine cnc_df_obs

  ! -- Procedure related to time series

  !> @brief Procedure related to time series
  !!
  !! Assign tsLink%Text appropriately for all time series in use by package.
  !! For any specified dependent variable package, for example either the
  !! constant concentration or constant temperature packages, the dependent
  !! variable can be controlled by time series.
  !<
  subroutine cnc_rp_ts(this)
    ! -- dummy
    class(GwtCncType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, nlinks
    type(TimeSeriesLinkType), pointer :: tslink => null()
    !
    nlinks = this%TsManager%boundtslinks%Count()
    do i = 1, nlinks
      tslink => GetTimeSeriesLinkFromList(this%TsManager%boundtslinks, i)
      if (associated(tslink)) then
        select case (tslink%JCol)
        case (1)
          tslink%Text = trim(this%depvartype)
        end select
      end if
    end do
  end subroutine cnc_rp_ts

  !> @brief Apply auxiliary multiplier to specified concentration if
  !< appropriate
  function conc_mult(this, row) result(conc)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(GwtCncType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: conc
    !
    if (this%iauxmultcol > 0) then
      conc = this%tspvar(row) * this%auxvar(this%iauxmultcol, row)
    else
      conc = this%tspvar(row)
    end if
  end function conc_mult

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index and row.
  !<
  function cnc_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(GwtCncType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%conc_mult(row)
    case default
      write (errmsg, '(3a)') 'Programming error. ', &
               & adjustl(trim(this%filtyp)), ' bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function cnc_bound_value

end module GwtCncModule
