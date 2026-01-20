module DrnModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, DTWO, &
                             LENFTYPE, LENPACKAGENAME, LENAUXNAME, LINELENGTH
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryHelperModule, only: create_mem_path
  use SmoothingModule, only: sQSaturation, sQSaturationDerivative, &
                             sQuadraticSaturation
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: drn_create
  public :: DrnType
  !
  character(len=LENFTYPE) :: ftype = 'DRN'
  character(len=LENPACKAGENAME) :: text = '             DRN'
  !
  type, extends(BndExtType) :: DrnType

    real(DP), dimension(:), pointer, contiguous :: elev => null() !< DRN elevation
    real(DP), dimension(:), pointer, contiguous :: cond => null() !< DRN conductance at aquifer interface
    integer(I4B), pointer :: iauxddrncol => null()
    integer(I4B), pointer :: icubic_scaling => null()

  contains

    procedure :: allocate_scalars => drn_allocate_scalars
    procedure :: allocate_arrays => drn_allocate_arrays
    procedure :: source_options => drn_options
    procedure :: log_drn_options
    procedure :: bnd_rp => drn_rp
    procedure :: bnd_ck => drn_ck
    procedure :: bnd_cf => drn_cf
    procedure :: bnd_fc => drn_fc
    procedure :: bnd_fn => drn_fn
    procedure :: bnd_da => drn_da
    procedure :: define_listlabel
    procedure :: get_drain_elevations
    procedure :: get_drain_factor
    procedure :: bound_value => drn_bound_value
    procedure :: cond_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => drn_obs_supported
    procedure, public :: bnd_df_obs => drn_df_obs
    procedure, public :: drn_store_user_cond
  end type DrnType

contains

  !> @brief Create a New Drn Package and point packobj to the new package
  !<
  subroutine drn_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(DrnType), pointer :: drnobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (drnobj)
    packobj => drnobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call drnobj%allocate_scalars()
    !s
    ! -- initialize package
    call packobj%pack_initialize()
    !
    ! -- initialize
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine drn_create

  !> @brief Deallocate memory
  !<
  subroutine drn_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(DrnType) :: this
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- scalars
    call mem_deallocate(this%iauxddrncol)
    call mem_deallocate(this%icubic_scaling)
    !
    ! -- arrays
    call mem_deallocate(this%elev, 'ELEV', this%memoryPath)
    call mem_deallocate(this%cond, 'COND', this%memoryPath)
  end subroutine drn_da

  !> @brief Allocate package scalar members
  !<
  subroutine drn_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(DrnType) :: this
    !
    ! -- call base type allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%iauxddrncol, 'IAUXDDRNCOL', this%memoryPath)
    call mem_allocate(this%icubic_scaling, 'ICUBIC_SCALING', this%memoryPath)
    !
    ! -- Set values
    this%iauxddrncol = 0
    if (this%inewton /= 0) then
      this%icubic_scaling = 1
    else
      this%icubic_scaling = 0
    end if
  end subroutine drn_allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine drn_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(DrnType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    !
    ! -- call base type allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set drn input context pointers
    call mem_setptr(this%elev, 'ELEV', this%input_mempath)
    call mem_setptr(this%cond, 'COND', this%input_mempath)
    !
    ! --checkin drn input context pointers
    call mem_checkin(this%elev, 'ELEV', this%memoryPath, &
                     'ELEV', this%input_mempath)
    call mem_checkin(this%cond, 'COND', this%memoryPath, &
                     'COND', this%input_mempath)
  end subroutine drn_allocate_arrays

  !> @brief Read and prepare
  !<
  subroutine drn_rp(this)
    use TdisModule, only: kper
    ! -- dummy
    class(DrnType), intent(inout) :: this
    !
    if (this%iper /= kper) return
    !
    ! -- Call the parent class read and prepare
    call this%BndExtType%bnd_rp()
    !
    ! -- store user cond
    if (this%ivsc == 1) then
      call this%drn_store_user_cond()
    end if
    !
    ! -- Write the list to iout if requested
    if (this%iprpak /= 0) then
      call this%write_list()
    end if
  end subroutine drn_rp

  !> @brief Source options specific to DrnType
  !<
  subroutine drn_options(this)
    ! -- modules
    use InputOutputModule, only: urword
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use GwfDrnInputModule, only: GwfDrnParamFoundType
    ! -- dummy
    class(DrnType), intent(inout) :: this
    ! -- local
    type(GwfDrnParamFoundType) :: found
    character(len=LENAUXNAME) :: ddrnauxname
    integer(I4B) :: n
    !
    ! -- source base class options
    call this%BndExtType%source_options()
    !
    ! -- source drain options
    call mem_set_value(this%imover, 'MOVER', this%input_mempath, found%mover)
    call mem_set_value(ddrnauxname, 'AUXDEPTHNAME', this%input_mempath, &
                       found%auxdepthname)
    call mem_set_value(this%icubic_scaling, 'ICUBICSFAC', this%input_mempath, &
                       found%icubicsfac)
    !
    if (found%auxdepthname) then
      this%iauxddrncol = -1
      !
      ! -- Error if no aux variable specified
      if (this%naux == 0) then
        write (errmsg, '(a,2(1x,a))') &
          'AUXDEPTHNAME was specified as', trim(adjustl(ddrnauxname)), &
          'but no AUX variables specified.'
        call store_error(errmsg)
      end if
      !
      ! -- Assign ddrn column
      this%iauxddrncol = 0
      do n = 1, this%naux
        if (ddrnauxname == this%auxname(n)) then
          this%iauxddrncol = n
          exit
        end if
      end do
      !
      ! -- Error if aux variable cannot be found
      if (this%iauxddrncol == 0) then
        write (errmsg, '(a,2(1x,a))') &
          'AUXDEPTHNAME was specified as', trim(adjustl(ddrnauxname)), &
          'but no AUX variable found with this name.'
        call store_error(errmsg)
      end if
    end if
    !
    ! -- log DRN specific options
    call this%log_drn_options(found)
  end subroutine drn_options

  !> @ brief Log DRN specific package options
  !<
  subroutine log_drn_options(this, found)
    ! -- modules
    use GwfDrnInputModule, only: GwfDrnParamFoundType
    ! -- dummy variables
    class(DrnType), intent(inout) :: this !< BndExtType object
    type(GwfDrnParamFoundType), intent(in) :: found
    ! -- local variables
    ! -- format
    !
    ! -- log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'
    !
    if (found%mover) then
      write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
    end if
    !
    if (found%icubicsfac) then
      write (this%iout, '(4x,a,1x,a)') &
        'CUBIC SCALING will be used for drains with non-zero DDRN values', &
        'even if the NEWTON-RAPHSON method is not being used.'
    end if
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine log_drn_options

  !> @brief Check drain boundary condition data
  !<
  subroutine drn_ck(this)
    ! -- dummy
    class(DrnType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: node
    real(DP) :: bt
    real(DP) :: drndepth
    real(DP) :: drntop
    real(DP) :: drnbot
    ! -- formats
    character(len=*), parameter :: fmtddrnerr = &
      "('SCALED-CONDUCTANCE DRN BOUNDARY (',i0,') BOTTOM ELEVATION &
      &(',f10.3,') IS LESS THAN CELL BOTTOM (',f10.3,')')"
    character(len=*), parameter :: fmtdrnerr = &
      "('DRN BOUNDARY (',i0,') ELEVATION (',f10.3,') IS LESS THAN CELL &
      &BOTTOM (',f10.3,')')"
    character(len=*), parameter :: fmtcondmulterr = &
      "('DRN BOUNDARY (',i0,') CONDUCTANCE MULTIPLIER (',g10.3,') IS &
      &LESS THAN ZERO')"
    character(len=*), parameter :: fmtconderr = &
      "('DRN BOUNDARY (',i0,') CONDUCTANCE (',g10.3,') IS LESS THAN &
      &ZERO')"
    !
    ! -- check stress period data
    do i = 1, this%nbound
      node = this%nodelist(i)
      bt = this%dis%bot(node)
      !
      ! -- calculate the drainage depth and the top and bottom of
      !    the conductance scaling elevations
      call this%get_drain_elevations(i, drndepth, drntop, drnbot)
      !
      ! -- accumulate errors
      if (drnbot < bt .and. this%icelltype(node) /= 0) then
        if (drndepth /= DZERO) then
          write (errmsg, fmt=fmtddrnerr) i, drnbot, bt
        else
          write (errmsg, fmt=fmtdrnerr) i, drnbot, bt
        end if
        call store_error(errmsg)
      end if
      if (this%iauxmultcol > 0) then
        if (this%auxvar(this%iauxmultcol, i) < DZERO) then
          write (errmsg, fmt=fmtcondmulterr) &
            i, this%auxvar(this%iauxmultcol, i)
          call store_error(errmsg)
        end if
      end if
      if (this%cond(i) < DZERO) then
        write (errmsg, fmt=fmtconderr) i, this%cond(i)
        call store_error(errmsg)
      end if
    end do
    !
    ! -- write summary of drain package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine drn_ck

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Skip if no drains
  !<
  subroutine drn_cf(this)
    ! -- dummy
    class(DrnType) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: node
    real(DP) :: cdrn
    real(DP) :: drnbot
    real(DP) :: fact
    !
    ! -- Return if no drains
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each drn entry
    do i = 1, this%nbound
      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if
      !
      ! -- set local variables for this drain
      cdrn = this%cond_mult(i)

      !
      ! -- calculate the drainage scaling factor
      call this%get_drain_factor(i, fact, drnbot)
      !
      ! -- calculate rhs and hcof
      this%rhs(i) = -fact * cdrn * drnbot
      this%hcof(i) = -fact * cdrn
    end do
  end subroutine drn_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine drn_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(DrnType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: ipos
    real(DP) :: fact
    real(DP) :: drnbot
    real(DP) :: drncond
    real(DP) :: qdrn
    !
    ! -- packmvrobj fc
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
      ! -- calculate the drainage scaling factor
      call this%get_drain_factor(i, fact, drnbot)
      !
      ! -- If mover is active and this drain is discharging,
      !    store available water (as positive value).
      if (this%imover == 1 .and. fact > DZERO) then
        drncond = this%cond_mult(i)
        qdrn = fact * drncond * (this%xnew(n) - drnbot)
        call this%pakmvrobj%accumulate_qformvr(i, qdrn)
      end if
    end do
  end subroutine drn_fc

  !> @brief Fill newton terms
  !<
  subroutine drn_fn(this, rhs, ia, idxglo, matrix_sln)
    implicit none
    ! -- dummy
    class(DrnType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: ipos
    real(DP) :: cdrn
    real(DP) :: xnew
    real(DP) :: drndepth
    real(DP) :: drntop
    real(DP) :: drnbot
    real(DP) :: drterm
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    if (this%iauxddrncol /= 0) then
      do i = 1, this%nbound
        node = this%nodelist(i)
        !
        ! -- test if node is constant or inactive
        if (this%ibound(node) <= 0) then
          cycle
        end if
        !
        ! -- set local variables for this drain
        cdrn = this%cond_mult(i)
        xnew = this%xnew(node)
        !
        ! -- calculate the drainage depth and the top and bottom of
        !    the conductance scaling elevations
        call this%get_drain_elevations(i, drndepth, drntop, drnbot)
        !
        ! -- calculate scaling factor
        if (drndepth /= DZERO) then
          drterm = sQSaturationDerivative(drntop, drnbot, xnew, &
                                          c1=-DONE, c2=DTWO)
          drterm = drterm * cdrn * (drnbot - xnew)
          !
          ! -- fill amat and rhs with newton-raphson terms
          ipos = ia(node)
          call matrix_sln%add_value_pos(idxglo(ipos), drterm)
          rhs(node) = rhs(node) + drterm * xnew
        end if
      end do
    end if
  end subroutine drn_fn

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(DrnType), intent(inout) :: this
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'DRAIN EL.'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'CONDUCTANCE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  !> @brief Define drain depth and the top and bottom elevations used to scale
  !! the drain conductance
  !<
  subroutine get_drain_elevations(this, i, drndepth, drntop, drnbot)
    ! -- dummy
    class(DrnType), intent(inout) :: this
    integer(I4B), intent(in) :: i
    real(DP), intent(inout) :: drndepth
    real(DP), intent(inout) :: drntop
    real(DP), intent(inout) :: drnbot
    ! -- local
    real(DP) :: drnelev
    real(DP) :: elev
    !
    ! -- initialize dummy and local variables
    drndepth = DZERO
    drnelev = this%elev(i)
    !
    ! -- set the drain depth
    if (this%iauxddrncol > 0) then
      drndepth = this%auxvar(this%iauxddrncol, i)
    end if
    !
    ! -- calculate the top and bottom drain elevations
    if (drndepth /= DZERO) then
      elev = drnelev + drndepth
      drntop = max(elev, drnelev)
      drnbot = min(elev, drnelev)
    else
      drntop = drnelev
      drnbot = drnelev
    end if
  end subroutine get_drain_elevations

  !> @brief Get the drain conductance scale factor
  !<
  subroutine get_drain_factor(this, i, factor, opt_drnbot)
    ! -- dummy
    class(DrnType), intent(inout) :: this
    integer(I4B), intent(in) :: i
    real(DP), intent(inout) :: factor
    real(DP), intent(inout), optional :: opt_drnbot
    ! -- local
    integer(I4B) :: node
    real(DP) :: xnew
    real(DP) :: drndepth
    real(DP) :: drntop
    real(DP) :: drnbot
    !
    ! -- set local variables for this drain
    node = this%nodelist(i)
    xnew = this%xnew(node)
    !
    ! -- calculate the drainage depth and the top and bottom of
    !    the conductance scaling elevations
    call this%get_drain_elevations(i, drndepth, drntop, drnbot)
    !
    ! -- set opt_drnbot to drnbot if passed as dummy variable
    if (present(opt_drnbot)) then
      opt_drnbot = drnbot
    end if
    !
    ! -- calculate scaling factor
    if (drndepth /= DZERO) then
      if (this%icubic_scaling /= 0) then
        factor = sQSaturation(drntop, drnbot, xnew, c1=-DONE, c2=DTWO)
      else
        factor = sQuadraticSaturation(drntop, drnbot, xnew, eps=DZERO)
      end if
    else
      if (xnew <= drnbot) then
        factor = DZERO
      else
        factor = DONE
      end if
    end if
  end subroutine get_drain_factor

  ! -- Procedures related to observations

  !> @brief Return true because DRN package supports observations
  !!
  !! Overrides BndType%bnd_obs_supported()
  !<
  logical function drn_obs_supported(this)
    implicit none
    ! -- dummy
    class(DrnType) :: this
    !
    drn_obs_supported = .true.
  end function drn_obs_supported

  !> @brief Store observation type supported by DRN package
  !!
  !! Overrides BndType%bnd_df_obs
  !<
  subroutine drn_df_obs(this)
    implicit none
    ! -- dummy
    class(DrnType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType('drn', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine drn_df_obs

  !> @brief Store user-specified drain conductance
  !<
  subroutine drn_store_user_cond(this)
    ! -- dummy
    class(DrnType), intent(inout) :: this !< BndExtType object
    ! -- local
    integer(I4B) :: n
    !
    ! -- store backup copy of conductance values
    do n = 1, this%nbound
      this%condinput(n) = this%cond_mult(n)
    end do
  end subroutine drn_store_user_cond

  !> @brief Apply multiplier to conductance value depending on user-selected option
  !<
  function cond_mult(this, row) result(cond)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(DrnType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: cond
    !
    if (this%iauxmultcol > 0) then
      cond = this%cond(row) * this%auxvar(this%iauxmultcol, row)
    else
      cond = this%cond(row)
    end if
  end function cond_mult

  !> @brief Return requested boundary value
  !<
  function drn_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(DrnType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%elev(row)
    case (2)
      bndval = this%cond_mult(row)
    case default
      errmsg = 'Programming error. DRN bound value requested column '&
               &'outside range of ncolbnd (2).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function drn_bound_value

end module DrnModule
