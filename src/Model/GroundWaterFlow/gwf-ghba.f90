module ghbamodule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DNODATA, LENFTYPE, LENPACKAGENAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: ghba_create
  public :: GhbaType
  !
  character(len=LENFTYPE) :: ftype = 'GHBA'
  character(len=LENPACKAGENAME) :: text = '            GHBA'
  !
  type, extends(BndExtType) :: GhbaType
    real(DP), dimension(:), pointer, contiguous :: bhead => null() !< GHBA boundary head
    real(DP), dimension(:), pointer, contiguous :: cond => null() !< GHBA hydraulic conductance
  contains
    procedure :: allocate_arrays => ghba_allocate_arrays
    procedure :: source_options => ghba_options
    procedure :: source_dimensions => ghba_dimensions
    procedure :: log_ghba_options
    procedure :: bnd_rp => ghba_rp
    procedure :: bnd_ck => ghba_ck
    procedure :: bnd_cf => ghba_cf
    procedure :: bnd_fc => ghba_fc
    procedure :: bnd_da => ghba_da
    procedure :: define_listlabel
    procedure :: bound_value => ghba_bound_value
    procedure :: cond_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => ghba_obs_supported
    procedure, public :: bnd_df_obs => ghba_df_obs
    procedure, public :: ghba_store_user_cond
  end type GhbaType

contains

  !> @brief Create a New Ghba Package and point bndobj to the new package
  !<
  subroutine ghba_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(GhbaType), pointer :: ghbaobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (ghbaobj)
    packobj => ghbaobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call packobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine ghba_create

  !> @brief Deallocate memory
  !<
  subroutine ghba_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GhbaType) :: this
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- arrays
    call mem_deallocate(this%bhead, 'BHEAD', this%memoryPath)
    call mem_deallocate(this%cond, 'COND', this%memoryPath)
  end subroutine ghba_da

  !> @brief Set options specific to GhbaType
  !<
  subroutine ghba_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use GwfGhbaInputModule, only: GwfGhbaParamFoundType
    ! -- dummy
    class(GhbaType), intent(inout) :: this
    ! -- local
    type(GwfGhbaParamFoundType) :: found
    !
    ! -- source base class options
    call this%BndExtType%source_options()
    !
    ! -- source options from input context
    call mem_set_value(this%imover, 'MOVER', this%input_mempath, found%mover)
    !
    ! -- log ghba specific options
    call this%log_ghba_options(found)
  end subroutine ghba_options

  !> @brief Log options specific to GhbaType
  !<
  subroutine log_ghba_options(this, found)
    ! -- modules
    use GwfGhbaInputModule, only: GwfGhbaParamFoundType
    ! -- dummy
    class(GhbaType), intent(inout) :: this !< BndExtType object
    type(GwfGhbaParamFoundType), intent(in) :: found
    !
    ! -- log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'
    !
    if (found%mover) then
      write (this%iout, '(4x,A)') 'MOVER OPTION ENABLED'
    end if
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine log_ghba_options

  !> @brief Set dimensions specific to GhbaType
  !<
  subroutine ghba_dimensions(this)
    ! -- modules
    ! -- dummy
    class(GhbaType), intent(inout) :: this
    ! -- local
    !
    ! -- set maxbound
    this%maxbound = this%dis%nodesuser
    !
    ! -- set nbound, which applies for duration of simulation
    this%nbound = this%dis%nodesuser
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
  end subroutine ghba_dimensions

  !> @brief Allocate arrays
  !<
  subroutine ghba_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(GhbaType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    !
    ! -- call base type allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set ghba input context pointers
    call mem_setptr(this%bhead, 'BHEAD', this%input_mempath)
    call mem_setptr(this%cond, 'COND', this%input_mempath)
    !
    ! --checkin ghba input context pointers
    call mem_checkin(this%bhead, 'BHEAD', this%memoryPath, &
                     'BHEAD', this%input_mempath)
    call mem_checkin(this%cond, 'COND', this%memoryPath, &
                     'COND', this%input_mempath)
  end subroutine ghba_allocate_arrays

  !> @brief Read and prepare
  !<
  subroutine ghba_rp(this)
    ! -- modules
    use TdisModule, only: kper
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(GhbaType), intent(inout) :: this
    integer(I4B) :: i, noder
    character(len=LINELENGTH) :: nodestr
    !
    if (this%iper /= kper) return
    !
    ! -- Update the nodelist
    do i = 1, this%nbound
      if (this%bhead(i) == DNODATA) then
        this%nodelist(i) = 0
      else
        noder = this%dis%get_nodenumber(i, 1)
        if (noder > 0) then
          this%nodelist(i) = noder
        else
          call this%dis%nodeu_to_string(i, nodestr)
          write (errmsg, *) &
            ' Cell is outside active grid domain: '// &
            trim(adjustl(nodestr))
          call store_error(errmsg)
        end if
      end if
    end do
    !
    ! -- exit if errors were found
    if (count_errors() > 0) then
      write (errmsg, *) count_errors(), ' errors encountered.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- store user cond
    if (this%ivsc == 1) then
      call this%ghba_store_user_cond()
    end if
    !
    ! -- Write the list to iout if requested
    if (this%iprpak /= 0) then
      call this%write_list()
    end if
  end subroutine ghba_rp

  !> @brief Check ghba boundary condition data
  !<
  subroutine ghba_ck(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error, count_errors, store_error_unit
    ! -- dummy
    class(GhbaType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: i
    integer(I4B) :: noder
    real(DP) :: bt
    ! -- formats
    character(len=*), parameter :: fmtghberr = &
      "('GHBA BOUNDARY (',i0,') HEAD (',f10.3,') IS LESS THAN CELL &
      &BOTTOM (',f10.3,')')"
    character(len=*), parameter :: fmtcondmulterr = &
      "('GHBA BOUNDARY (',i0,') CONDUCTANCE MULTIPLIER (',g10.3,') IS &
      &NO DATA VALUE OR LESS THAN ZERO')"
    character(len=*), parameter :: fmtconderr = &
      "('GHBA BOUNDARY (',i0,') CONDUCTANCE (',g10.3,') IS NO DATA VALUE &
      &OR LESS THAN ZERO')"
    !
    ! -- check stress period data
    do i = 1, this%nbound
      noder = this%nodelist(i)
      if (noder == 0) cycle
      bt = this%dis%bot(noder)
      ! -- accumulate errors
      if (this%bhead(i) < bt .and. this%icelltype(noder) /= 0) then
        write (errmsg, fmt=fmtghberr) i, this%bhead(i), bt
        call store_error(errmsg)
      end if
      if (this%iauxmultcol > 0) then
        if (this%auxvar(this%iauxmultcol, i) == DNODATA .or. &
            this%auxvar(this%iauxmultcol, i) < DZERO) then
          write (errmsg, fmt=fmtcondmulterr) &
            i, this%auxvar(this%iauxmultcol, i)
          call store_error(errmsg)
        end if
      end if
      if (this%cond(i) == DNODATA .or. &
          this%cond(i) < DZERO) then
        write (errmsg, fmt=fmtconderr) i, this%cond(i)
        call store_error(errmsg)
      end if
    end do
    !
    !write summary of ghb package error messages
    if (count_errors() > 0) then
      call store_error_unit(this%inunit)
    end if
  end subroutine ghba_ck

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! Skip if no GHBs
  !<
  subroutine ghba_cf(this)
    ! -- dummy
    class(GhbaType) :: this
    ! -- local
    integer(I4B) :: i, noder
    !
    do i = 1, this%nbound
      noder = this%nodelist(i)
      if (noder == 0) cycle
      if (this%ibound(noder) .le. 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if
      this%hcof(i) = -this%cond_mult(i)
      this%rhs(i) = -this%cond_mult(i) * this%bhead(i)
    end do
  end subroutine ghba_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine ghba_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GhbaType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: i, noder, ipos
    real(DP) :: cond, bhead, qghba
    !
    ! -- pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if

    do i = 1, this%nbound
      noder = this%nodelist(i)
      if (noder == 0) cycle
      rhs(noder) = rhs(noder) + this%rhs(i)
      ipos = ia(noder)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
      !
      ! -- If mover is active and this boundary is discharging,
      !    store available water (as positive value).
      bhead = this%bhead(i)
      if (this%imover == 1 .and. this%xnew(noder) > bhead) then
        cond = this%cond_mult(i)
        qghba = cond * (this%xnew(noder) - bhead)
        call this%pakmvrobj%accumulate_qformvr(i, qghba)
      end if
    end do
  end subroutine ghba_fc

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(GhbaType), intent(inout) :: this
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'STAGE'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'CONDUCTANCE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  ! -- Procedures related to observations

  !> @brief Return true because GHBA package supports observations
  !!
  !! Overrides BndType%bnd_obs_supported()
  !<
  logical function ghba_obs_supported(this)
    implicit none
    ! -- dummy
    class(GhbaType) :: this
    !
    ghba_obs_supported = .true.
  end function ghba_obs_supported

  !> @brief Store observation type supported by GHBA package
  !!
  !! Overrides BndType%bnd_df_obs
  !<
  subroutine ghba_df_obs(this)
    implicit none
    ! -- dummy
    class(GhbaType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType('ghba', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine ghba_df_obs

  !> @brief Store user-specified conductance for GHBA boundary type
  !<
  subroutine ghba_store_user_cond(this)
    ! -- modules
    ! -- dummy
    class(GhbaType), intent(inout) :: this !< BndExtType object
    ! -- local
    integer(I4B) :: n
    !
    ! -- store backup copy of conductance values
    do n = 1, this%nbound
      this%condinput(n) = this%cond_mult(n)
    end do
  end subroutine ghba_store_user_cond

  !> @brief Apply multiplier to GHBA conductance if option AUXMULTCOL is used
  !<
  function cond_mult(this, row) result(cond)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(GhbaType), intent(inout) :: this !< BndExtType object
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
  function ghba_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(GhbaType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%bhead(row)
    case (2)
      bndval = this%cond_mult(row)
    case default
      errmsg = 'Programming error. GHBA bound value requested column '&
               &'outside range of ncolbnd (2).'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end select
  end function ghba_bound_value

end module ghbamodule
