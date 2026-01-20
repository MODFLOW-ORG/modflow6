module ChdModule
  !
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, NAMEDBOUNDFLAG, LENFTYPE, &
                             LINELENGTH, LENPACKAGENAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryHelperModule, only: create_mem_path
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
  public :: chd_create, ChdType
  !
  character(len=LENFTYPE) :: ftype = 'CHD'
  character(len=LENPACKAGENAME) :: text = '             CHD'
  !
  type, extends(BndExtType) :: ChdType
    real(DP), dimension(:), pointer, contiguous :: head => null() !< constant head array
    real(DP), dimension(:), pointer, contiguous :: ratechdin => null() !simulated flows into constant head (excluding other chds)
    real(DP), dimension(:), pointer, contiguous :: ratechdout => null() !simulated flows out of constant head (excluding to other chds)
  contains
    procedure :: bnd_rp => chd_rp
    procedure :: bnd_ad => chd_ad
    procedure :: bnd_ck => chd_ck
    procedure :: bnd_fc => chd_fc
    procedure :: bnd_cq => chd_cq
    procedure :: bnd_bd => chd_bd
    procedure :: bnd_da => chd_da
    procedure :: allocate_arrays => chd_allocate_arrays
    procedure :: define_listlabel
    procedure :: bound_value => chd_bound_value
    procedure :: head_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => chd_obs_supported
    procedure, public :: bnd_df_obs => chd_df_obs
    !
    procedure, private :: calc_chd_rate
  end type ChdType

contains

  !> @brief Create a new constant head package
  !!
  !! Routine points packobj to the newly created package
  !<
  subroutine chd_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    ! -- local
    type(ChdType), pointer :: chdobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (chdobj)
    packobj => chdobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call chdobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    ! -- store values
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine chd_create

  !> @brief Allocate arrays specific to the constant head package
  !<
  subroutine chd_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(ChdType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    integer(I4B) :: i
    !
    ! -- call standard BndType allocate scalars
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- allocate ratechdex
    call mem_allocate(this%ratechdin, this%maxbound, 'RATECHDIN', this%memoryPath)
    call mem_allocate(this%ratechdout, this%maxbound, 'RATECHDOUT', &
                      this%memoryPath)
    do i = 1, this%maxbound
      this%ratechdin(i) = DZERO
      this%ratechdout(i) = DZERO
    end do
    !
    ! -- set constant head array input context pointer
    call mem_setptr(this%head, 'HEAD', this%input_mempath)
    !
    ! -- checkin constant head array input context pointer
    call mem_checkin(this%head, 'HEAD', this%memoryPath, &
                     'HEAD', this%input_mempath)
  end subroutine chd_allocate_arrays

  !> @brief Constant concentration/temperature read and prepare (rp) routine
  !<
  subroutine chd_rp(this)
    !
    use TdisModule, only: kper
    ! -- dummy
    class(ChdType), intent(inout) :: this
    ! -- local
    character(len=30) :: nodestr
    integer(I4B) :: i, node, ibd, ierr
    !
    if (this%iper /= kper) return
    !
    ! -- Reset previous CHDs to active cell
    do i = 1, this%nbound
      node = this%nodelist(i)
      this%ibound(node) = this%ibcnum
    end do
    !
    ! -- Call the parent class read and prepare
    call this%BndExtType%bnd_rp()
    !
    ! -- Set ibound to -(ibcnum + 1) for constant head cells
    ierr = 0
    do i = 1, this%nbound
      node = this%nodelist(i)
      ibd = this%ibound(node)
      if (ibd < 0) then
        call this%dis%noder_to_string(node, nodestr)
        write (errmsg, '(3a)') &
          'Cell is already a constant head (', trim(adjustl(nodestr)), ').'
        call store_error(errmsg)
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
  end subroutine chd_rp

  !> @brief Constant head package advance routine
  !!
  !! Add package connections to matrix
  !<
  subroutine chd_ad(this)
    ! -- modules
    ! -- dummy
    class(ChdType) :: this
    ! -- local
    integer(I4B) :: i, node
    real(DP) :: hb
    ! -- formats
    !
    ! -- Process each entry in the specified-head cell list
    do i = 1, this%nbound
      node = this%nodelist(i)
      hb = this%head_mult(i)
      !
      this%xnew(node) = hb
      this%xold(node) = this%xnew(node)
    end do
    !
    ! -- For each observation, push simulated value and corresponding
    !    simulation time from "current" to "preceding" and reset
    !    "current" value.
    call this%obs%obs_ad()
  end subroutine chd_ad

  !> @brief Check constant concentration/temperature boundary condition data
  !<
  subroutine chd_ck(this)
    ! -- modules
    ! -- dummy
    class(ChdType), intent(inout) :: this
    ! -- local
    character(len=30) :: nodestr
    integer(I4B) :: i
    integer(I4B) :: node
    real(DP) :: bt
    ! -- formats
    character(len=*), parameter :: fmtchderr = &
      "('CHD BOUNDARY ',i0,' HEAD (',g0,') IS LESS THAN CELL &
      &BOTTOM (',g0,')',' FOR CELL ',a)"
    !
    ! -- check stress period data
    do i = 1, this%nbound
      node = this%nodelist(i)
      bt = this%dis%bot(node)
      ! -- accumulate errors
      if (this%head_mult(i) < bt .and. this%icelltype(node) /= 0) then
        call this%dis%noder_to_string(node, nodestr)
        write (errmsg, fmt=fmtchderr) i, this%head_mult(i), bt, trim(nodestr)
        call store_error(errmsg)
      end if
    end do
    !
    ! write summary of chd package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine chd_ck

  !> @brief Override bnd_fc and do nothing
  !!
  !! For constant head boundary type, the call to bnd_fc
  !! needs to be overwritten to do nothing
  !<
  subroutine chd_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(ChdType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
  end subroutine chd_fc

  !> @brief Calculate flow associated with constant head boundary
  !!
  !! This method overrides bnd_cq()
  !<
  subroutine chd_cq(this, x, flowja, iadv)
    class(ChdType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv

    ! NB: the rate calculation cannot be done until chd_bd below

  end subroutine chd_cq

  !> @brief Calculate the CHD cell rates, to be called
  !< after all updates to the model flowja are done
  subroutine calc_chd_rate(this)
    ! -- modules
    ! -- dummy
    class(ChdType), intent(inout) :: this
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
          q = this%flowja(ipos)
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
        ! -- For chd, store total flow in rhs so it is available for other
        !    calculations
        this%rhs(i) = -rate
        this%hcof(i) = DZERO
        !
        ! -- Save simulated value to simvals array.
        this%simvals(i) = rate
        this%ratechdin(i) = ratein
        this%ratechdout(i) = rateout
        this%flowja(idiag) = this%flowja(idiag) + rate
        !
      end do
      !
    end if
  end subroutine calc_chd_rate

  !> @brief Add package ratin/ratout to model budget
  !<
  subroutine chd_bd(this, model_budget)
    ! -- add package ratin/ratout to model budget
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    class(ChdType) :: this
    type(BudgetType), intent(inout) :: model_budget
    real(DP) :: ratin
    real(DP) :: ratout
    real(DP) :: dum
    integer(I4B) :: isuppress_output

    ! For CHDs at an exchange, under some conditions
    ! (XT3D), the model flowja into the cell is not
    ! finalized until after exg_cq. So we calculate
    ! the CHD rate here
    call this%calc_chd_rate()

    isuppress_output = 0
    call rate_accumulator(this%ratechdin(1:this%nbound), ratin, dum)
    call rate_accumulator(this%ratechdout(1:this%nbound), ratout, dum)
    call model_budget%addentry(ratin, ratout, delt, this%text, &
                               isuppress_output, this%packName)
  end subroutine chd_bd

  !> @brief Deallocate memory
  !<
  subroutine chd_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(ChdType) :: this
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- arrays
    call mem_deallocate(this%ratechdin)
    call mem_deallocate(this%ratechdout)
    call mem_deallocate(this%head, 'HEAD', this%memoryPath)
  end subroutine chd_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used.
  !<
  subroutine define_listlabel(this)
    class(ChdType), intent(inout) :: this
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'HEAD'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  ! -- Procedures related to observations

  !> @brief Overrides bnd_obs_supported from bndType class
  !!
  !! Return true since CHD package supports observations
  !<
  logical function chd_obs_supported(this)
    implicit none
    !
    class(ChdType) :: this
    !
    chd_obs_supported = .true.
  end function chd_obs_supported

  !> @brief Overrides bnd_df_obs from bndType class
  !!
  !! (1) Store observation type supported by CHD package and (2) override
  !! BndType%bnd_df_obs
  !<
  subroutine chd_df_obs(this)
    implicit none
    ! -- dummy
    class(ChdType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType('chd', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine chd_df_obs

  !> @brief Apply auxiliary multiplier to specified head if appropriate
  !<
  function head_mult(this, row) result(head)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(ChdType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: head
    !
    if (this%iauxmultcol > 0) then
      head = this%head(row) * this%auxvar(this%iauxmultcol, row)
    else
      head = this%head(row)
    end if
  end function head_mult

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index
  !!  and row.
  !<
  function chd_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(ChdType), intent(inout) :: this !< BndType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%head_mult(row)
    case default
      errmsg = 'Programming error. CHD bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function chd_bound_value

end module ChdModule
