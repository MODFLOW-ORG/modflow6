!> @brief This module contains the WEL package methods
!!
!! This module contains the overridden methods for the standard WEL package.
!! Several methods need to be overridden because of the AUTO_FLOW_REDUCE
!! option. Overridden methods include:
!!   - bnd_cf (AUTO_FLOW_REDUCE)
!!   - bnd_fc (AUTO_FLOW_REDUCE)
!!   - bnd_fn (AUTO_FLOW_REDUCE Newton-Raphson terms)
!!   - bnd_ot_package_flows (write AUTO_FLOW_REDUCE terms to csv file)
!!   - bnd_da (deallocate AUTO_FLOW_REDUCE variables)
!!   - bnd_bd_obs (wel-reduction observation added)
!!
!<

module WelModule
  ! -- modules used by WelModule methods
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DEM1, DONE, LENFTYPE, DNODATA, LINELENGTH
  use SimVariablesModule, only: errmsg, warnmsg
  use SimModule, only: store_error, store_error_filename, store_warning
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use SmoothingModule, only: sQSaturation, sQSaturationDerivative
  use ObserveModule, only: ObserveType
  use BlockParserModule, only: BlockParserType
  use InputOutputModule, only: GetUnit, openfile
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: wel_create
  !
  character(len=LENFTYPE) :: ftype = 'WEL' !< package ftype
  character(len=16) :: text = '             WEL' !< package flow text string
  !
  type, extends(BndExtType) :: WelType
    real(DP), dimension(:), pointer, contiguous :: q => null() !< volumetric well rate
    integer(I4B), pointer :: iflowred => null() !< flag indicating if the AUTO_FLOW_REDUCE option is active
    real(DP), pointer :: flowred => null() !< AUTO_FLOW_REDUCE variable
    integer(I4B), pointer :: ioutafrcsv => null() !< unit number for CSV output file containing wells with reduced puping rates
    integer(I4B), pointer :: iflowredlen => null() !< flag indicating flowred variable is a length value
  contains
    procedure :: allocate_scalars => wel_allocate_scalars
    procedure :: allocate_arrays => wel_allocate_arrays
    procedure :: source_options => wel_options
    procedure :: log_wel_options
    procedure :: bnd_rp => wel_rp
    procedure :: bnd_cf => wel_cf
    procedure :: bnd_fc => wel_fc
    procedure :: bnd_fn => wel_fn
    procedure :: bnd_da => wel_da
    procedure :: define_listlabel
    procedure :: bound_value => wel_bound_value
    procedure :: q_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => wel_obs_supported
    procedure, public :: bnd_df_obs => wel_df_obs
    procedure, public :: bnd_bd_obs => wel_bd_obs
    ! -- afr
    procedure, private :: wel_afr_csv_init
    procedure, private :: wel_afr_csv_write
  end type weltype

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new WEL Package object
  !!
  !<
  subroutine wel_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath)
    ! -- dummy variables
    class(BndType), pointer :: packobj !< pointer to default package type
    integer(I4B), intent(in) :: id !< package id
    integer(I4B), intent(in) :: ibcnum !< boundary condition number
    integer(I4B), intent(in) :: inunit !< unit number of WEL package input file
    integer(I4B), intent(in) :: iout !< unit number of model listing file
    character(len=*), intent(in) :: namemodel !< model name
    character(len=*), intent(in) :: pakname !< package name
    character(len=*), intent(in) :: mempath !< input mempath
    ! -- local variables
    type(WelType), pointer :: welobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (welobj)
    packobj => welobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call welobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine wel_create

  !> @ brief Deallocate package memory
    !!
    !!  Deallocate WEL package scalars and arrays.
    !!
  !<
  subroutine wel_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(WelType) :: this !< WelType object
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- scalars
    call mem_deallocate(this%iflowred)
    call mem_deallocate(this%flowred)
    call mem_deallocate(this%ioutafrcsv)
    call mem_deallocate(this%iflowredlen)
    call mem_deallocate(this%q, 'Q', this%memoryPath)
  end subroutine wel_da

  !> @ brief Allocate scalars
    !!
    !! Allocate and initialize scalars for the WEL package. The base model
    !! allocate scalars method is also called.
    !!
  !<
  subroutine wel_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(WelType) :: this !< WelType object
    !
    ! -- call base type allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%iflowred, 'IFLOWRED', this%memoryPath)
    call mem_allocate(this%flowred, 'FLOWRED', this%memoryPath)
    call mem_allocate(this%ioutafrcsv, 'IOUTAFRCSV', this%memoryPath)
    call mem_allocate(this%iflowredlen, 'IFLOWREDLEN', this%memoryPath)
    !
    ! -- Set values
    this%iflowred = 0
    this%ioutafrcsv = 0
    this%flowred = DZERO
    this%iflowredlen = 0
  end subroutine wel_allocate_scalars

  !> @ brief Allocate arrays
    !!
    !! Allocate and initialize arrays for the WEL package
    !!
  !<
  subroutine wel_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(WelType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    ! -- local
    !
    ! -- call BndExtType allocate scalars
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set constant head array input context pointer
    call mem_setptr(this%q, 'Q', this%input_mempath)
    !
    ! -- checkin constant head array input context pointer
    call mem_checkin(this%q, 'Q', this%memoryPath, &
                     'Q', this%input_mempath)
  end subroutine wel_allocate_arrays

  !> @ brief Source additional options for package
    !!
    !!  Source additional options for WEL package.
    !!
  !<
  subroutine wel_options(this)
    ! -- modules
    use InputOutputModule, only: urword
    use MemoryManagerExtModule, only: mem_set_value
    use GwfWelInputModule, only: GwfWelParamFoundType
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< WelType object
    ! -- local variables
    character(len=LINELENGTH) :: fname
    type(GwfWelParamFoundType) :: found
    ! -- formats
    character(len=*), parameter :: fmtflowred = &
      &"(4x, 'AUTOMATIC FLOW REDUCTION OF WELLS IMPLEMENTED.')"
    character(len=*), parameter :: fmtflowredv = &
      &"(4x, 'AUTOMATIC FLOW REDUCTION FRACTION (',g15.7,').')"
    !
    ! -- source base BndExtType options
    call this%BndExtType%source_options()
    !
    ! -- source well options from input context
    call mem_set_value(this%flowred, 'FLOWRED', this%input_mempath, found%flowred)
    call mem_set_value(fname, 'AFRCSVFILE', this%input_mempath, found%afrcsvfile)
    call mem_set_value(this%imover, 'MOVER', this%input_mempath, found%mover)
    call mem_set_value(this%iflowredlen, 'IFLOWREDLEN', this%input_mempath, &
                       found%iflowredlen)

    if (found%iflowredlen) then
      if (found%flowred .eqv. .FALSE.) then
        write (warnmsg, '(a)') &
          'FLOW_REDUCTION_LENGTH option specified but a AUTO_FLOW_REDUCTION value &
          &is not specified. The FLOW_REDUCTION_LENGTH option will be ignored.'
        call store_warning(warnmsg)
      else
        this%iflowredlen = 1
      end if
    end if

    if (found%flowred) then
      this%iflowred = 1
      if (this%flowred <= DZERO) then
        if (found%iflowredlen) then
          write (errmsg, '(a)') &
            'An AUTO_FLOW_REDUCTION value less than or equal to zero cannot be &
            &specified if the FLOW_REDUCTION_LENGTH option is specified.'
          call store_error(errmsg)
        else
          this%flowred = DEM1
        end if
      else if (this%flowred > DONE .and. this%iflowredlen == 0) then
        this%flowred = DONE
      end if
    end if

    if (found%afrcsvfile) then
      call this%wel_afr_csv_init(fname)
    end if

    if (found%mover) then
      this%imover = 1
    end if

    ! -- log WEL specific options
    call this%log_wel_options(found)
  end subroutine wel_options

  !> @ brief Log WEL specific package options
  !<
  subroutine log_wel_options(this, found)
    ! -- modules
    use GwfWelInputModule, only: GwfWelParamFoundType
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< BndExtType object
    type(GwfWelParamFoundType), intent(in) :: found
    ! -- local variables
    ! -- format
    character(len=*), parameter :: fmtflowred = &
      &"(4x, 'AUTOMATIC FLOW REDUCTION OF WELLS IMPLEMENTED.')"
    character(len=*), parameter :: fmtflowredv = &
      &"(4x, 'AUTOMATIC FLOW REDUCTION FRACTION (',g15.7,').')"
    character(len=*), parameter :: fmtflowredl = &
      &"(4x, 'AUTOMATIC FLOW REDUCTION LENGTH (',g15.7,').')"
    !
    ! -- log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'

    if (found%iflowredlen) then
      write (this%iout, fmtflowred)
      write (this%iout, '(4x,A)') &
        'AUTOMATIC FLOW REDUCTION FRACTION INTERPRETED AS A LENGTH'
    end if

    if (found%flowred) then
      if (this%iflowredlen == 0) then
        write (this%iout, fmtflowredv) this%flowred
      else
        write (this%iout, fmtflowredl) this%flowred
      end if
    end if
    !
    if (found%afrcsvfile) then
      ! -- currently no-op
    end if

    if (found%mover) then
      write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
    end if
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine log_wel_options

  !> @ brief WEL read and prepare
    !!
  !<
  subroutine wel_rp(this)
    use TdisModule, only: kper
    ! -- dummy
    class(WelType), intent(inout) :: this
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
  end subroutine wel_rp

  !> @ brief Formulate the package hcof and rhs terms.
    !!
    !!  Formulate the hcof and rhs terms for the WEL package that will be
    !!  added to the coefficient matrix and right-hand side vector.
    !!
  !<
  subroutine wel_cf(this)
    ! -- dummy variables
    class(WelType) :: this !< WelType object
    ! -- local variables
    integer(I4B) :: i, node, ict
    real(DP) :: qmult
    real(DP) :: q
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: thick
    !
    ! -- Return if no wells
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each well entry
    do i = 1, this%nbound
      node = this%nodelist(i)
      this%hcof(i) = DZERO
      if (this%ibound(node) <= 0) then
        this%rhs(i) = DZERO
        cycle
      end if
      q = this%q_mult(i)
      if (this%iflowred /= 0 .and. q < DZERO) then
        ict = this%icelltype(node)
        if (ict /= 0) then
          bt = this%dis%bot(node)
          if (this%iflowredlen == 0) then
            thick = this%dis%top(node) - bt
          else
            thick = DONE
          end if
          tp = bt + this%flowred * thick
          qmult = sQSaturation(tp, bt, this%xnew(node))
          q = q * qmult
        end if
      end if
      this%rhs(i) = -q
    end do
  end subroutine wel_cf

  !> @ brief Copy hcof and rhs terms into solution.
    !!
    !!  Add the hcof and rhs terms for the WEL package to the
    !!  coefficient matrix and right-hand side vector.
    !!
  !<
  subroutine wel_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy variables
    class(WelType) :: this !< WelType object
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
      ! -- If mover is active and this well is discharging,
      !    store available water (as positive value).
      if (this%imover == 1 .and. this%rhs(i) > DZERO) then
        call this%pakmvrobj%accumulate_qformvr(i, this%rhs(i))
      end if
    end do
  end subroutine wel_fc

  !> @ brief Add Newton-Raphson terms for package into solution.
    !!
    !!  Calculate and add the Newton-Raphson terms for the WEL package to the
    !!  coefficient matrix and right-hand side vector.
    !!
  !<
  subroutine wel_fn(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy variables
    class(WelType) :: this !< WelType object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! -- local variables
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: ipos
    integer(I4B) :: ict
    real(DP) :: drterm
    real(DP) :: q
    real(DP) :: tp
    real(DP) :: bt
    real(DP) :: thick
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      node = this%nodelist(i)
      !
      ! -- test if node is constant or inactive
      if (this%ibound(node) <= 0) then
        cycle
      end if
      !
      ! -- well rate is possibly head dependent
      ict = this%icelltype(node)
      if (this%iflowred /= 0 .and. ict /= 0) then
        ipos = ia(node)
        q = -this%rhs(i)
        if (q < DZERO) then
          ! -- calculate derivative for well
          tp = this%dis%top(node)
          bt = this%dis%bot(node)
          thick = tp - bt
          tp = bt + this%flowred * thick
          drterm = sQSaturationDerivative(tp, bt, this%xnew(node))
          drterm = drterm * this%q_mult(i)
          !--fill amat and rhs with newton-raphson terms
          call matrix_sln%add_value_pos(idxglo(ipos), drterm)
          rhs(node) = rhs(node) + drterm * this%xnew(node)
        end if
      end if
    end do
  end subroutine wel_fn

  !> @brief Initialize the auto flow reduce csv output file
  subroutine wel_afr_csv_init(this, fname)
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< WelType object
    character(len=*), intent(in) :: fname
    ! -- format
    character(len=*), parameter :: fmtafrcsv = &
      "(4x, 'AUTO FLOW REDUCE INFORMATION WILL BE SAVED TO FILE: ', a, /4x, &
      &'OPENED ON UNIT: ', I0)"

    this%ioutafrcsv = getunit()
    call openfile(this%ioutafrcsv, this%iout, fname, 'CSV', &
                  filstat_opt='REPLACE')
    write (this%iout, fmtafrcsv) trim(adjustl(fname)), &
      this%ioutafrcsv
    write (this%ioutafrcsv, '(a)') &
      'time,period,step,boundnumber,cellnumber,rate-requested,&
      &rate-actual,wel-reduction'
  end subroutine wel_afr_csv_init

  !> @brief Write out auto flow reductions only when & where they occur
  subroutine wel_afr_csv_write(this)
    ! -- modules
    use TdisModule, only: totim, kstp, kper
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< WelType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: nodereduced
    integer(I4B) :: nodeuser
    real(DP) :: v
    ! -- format
    do i = 1, this%nbound
      nodereduced = this%nodelist(i)
      !
      ! -- test if node is constant or inactive
      if (this%ibound(nodereduced) <= 0) then
        cycle
      end if
      v = this%q_mult(i) + this%rhs(i)
      if (v < DZERO) then
        nodeuser = this%dis%get_nodeuser(nodereduced)
        write (this%ioutafrcsv, '(*(G0,:,","))') &
          totim, kper, kstp, i, nodeuser, this%q_mult(i), this%simvals(i), v
      end if
    end do
  end subroutine wel_afr_csv_write

  !> @ brief Define the list label for the package
    !!
    !!  Method defined the list label for the WEL package. The list label is
    !!  the heading that is written to iout when PRINT_INPUT option is used.
    !!
  !<
  subroutine define_listlabel(this)
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< WelType object
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'STRESS RATE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  ! -- Procedures related to observations

  !> @brief Determine if observations are supported.
    !!
    !! Function to determine if observations are supported by the WEL package.
    !! Observations are supported by the WEL package.
    !!
    !! @return  wel_obs_supported       boolean indicating if observations are supported
    !!
  !<
  logical function wel_obs_supported(this)
    ! -- dummy variables
    class(WelType) :: this !< WelType object
    !
    ! -- set boolean
    wel_obs_supported = .true.
  end function wel_obs_supported

  !> @brief Define the observation types available in the package
    !!
    !! Method to define the observation types available in the WEL package.
    !!
  !<
  subroutine wel_df_obs(this)
    ! -- dummy variables
    class(WelType) :: this !< WelType object
    ! -- local variables
    integer(I4B) :: indx
    !
    ! -- initialize observations
    call this%obs%StoreObsType('wel', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for wel-reduction observation type.
    call this%obs%StoreObsType('wel-reduction', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine wel_df_obs

  !> @brief Save observations for the package
    !!
    !! Method to save simulated values for the WEL package.
    !!
  !<
  subroutine wel_bd_obs(this)
    ! -- dummy variables
    class(WelType) :: this !< WelType object
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
          case ('WEL')
            v = this%simvals(jj)
          case ('WEL-REDUCTION')
            if (this%iflowred > 0) then
              v = this%q_mult(jj) + this%rhs(jj)
            end if
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
    !
    ! -- Write the auto flow reduce csv file entries for this step
    if (this%ioutafrcsv > 0) then
      call this%wel_afr_csv_write()
    end if
  end subroutine wel_bd_obs

  function q_mult(this, row) result(q)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< BndExtType object
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
  function wel_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(WelType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%q_mult(row)
    case default
      errmsg = 'Programming error. WEL bound value requested column '&
               &'outside range of ncolbnd (1).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function wel_bound_value

end module WelModule
