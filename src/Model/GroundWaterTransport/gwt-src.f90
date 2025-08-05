module GwtSrcModule
  !
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DEM1, DONE, LENFTYPE, LENVARNAME
  use TspFmiModule, only: TspFmiType
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: src_create
  !
  character(len=LENFTYPE) :: ftype = 'SRC'
  character(len=16) :: text = '             SRC'
  !
  type, extends(BndExtType) :: GwtSrcType

    character(len=LENVARNAME) :: depvartype = '' !< stores string of dependent variable type, depending on model type
    type(TspFmiType), pointer :: fmi => null() ! pointer to GWE fmi object
    logical(LGP), pointer :: highest_sat => NULL()
    real(DP), dimension(:), pointer, contiguous :: smassrate => null() !< mass source loading rate
    integer(I4B), dimension(:), pointer, contiguous :: nodesontop => NULL() ! User provided cell numbers; nodelist is cells where recharge is applied)

  contains

    procedure :: allocate_scalars => src_allocate_scalars
    procedure :: allocate_arrays => src_allocate_arrays
    procedure :: source_options => src_options
    procedure :: bnd_rp => src_rp
    procedure :: bnd_cf => src_cf
    procedure :: bnd_fc => src_fc
    procedure :: bnd_da => src_da
    procedure :: define_listlabel
    procedure :: set_nodesontop
    procedure :: bound_value => src_bound_value
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => src_obs_supported
    procedure, public :: bnd_df_obs => src_df_obs

  end type GwtSrcType

contains

  !> @brief Create a source loading package
  !!
  !! This subroutine points bndobj to the newly created package
  !<
  subroutine src_create(packobj, id, ibcnum, inunit, iout, namemodel, &
                        depvartype, pakname, input_mempath, fmi)
    ! -- modules
    use BndModule, only: BndType
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=LENVARNAME), intent(in) :: depvartype
    character(len=*), intent(in) :: input_mempath
    type(TspFmiType), intent(in), target :: fmi
    ! -- local
    type(GwtSrcType), pointer :: srcobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (srcobj)
    packobj => srcobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, input_mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call srcobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1
    !
    ! -- Store the appropriate label based on the dependent variable
    srcobj%depvartype = depvartype

    srcobj%fmi => fmi

  end subroutine src_create

  !> @brief Set additional options specific to the GwtSrcType
  !<
  subroutine src_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use GwtSrcInputModule, only: GwtSrcParamFoundType
    ! -- dummy
    class(GwtSrcType), intent(inout) :: this
    ! -- local
    type(GwtSrcParamFoundType) :: found
    !
    ! -- source base class options
    call this%BndExtType%source_options()
    !
    ! -- source options from input context
    call mem_set_value(this%highest_sat, 'HIGHEST_SAT', this%input_mempath, &
                       found%highest_sat)

    write (this%iout, '(/1x,a)') 'PROCESSING SRC OPTIONS'
    if (found%highest_sat) then
      write (this%iout, '(4x,a)') &
        'Mass source loading rate will be applied to the highest cell at or below &
        &the specified cellid with a non-zero saturation.'
    end if
    write (this%iout, '(1x,a)') 'END OF SRC OPTIONS'
  end subroutine src_options

  !> @brief Deallocate memory
  !<
  subroutine src_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwtSrcType) :: this
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
    !
    ! -- arrays
    call mem_deallocate(this%smassrate, 'SMASSRATE', this%memoryPath)
    if (this%highest_sat) then
      call mem_deallocate(this%nodesontop, "NODESONTOP", this%memoryPath)
    end if
    !
    ! -- scalars
    call mem_deallocate(this%highest_sat)
  end subroutine src_da

  !> @brief Allocate scalars
  !!
  !! Allocate scalars specific to this source loading package
  !<
  subroutine src_allocate_scalars(this)
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtSrcType) :: this
    !
    ! -- base class allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%highest_sat, 'HIGHEST_SAT', this%memoryPath)
    !
    ! -- Set values
    this%highest_sat = .FALSE.

  end subroutine src_allocate_scalars

  !> @brief Allocate arrays
  !!
  !! Allocate scalars specific to this source loading package
  !<
  subroutine src_allocate_arrays(this, nodelist, auxvar)
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(GwtSrcType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist !< package nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar !< package aux variable array
    ! local
    integer(I4B) :: n
    !
    ! -- base class allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- allocate the object and assign values to object variables
    if (this%highest_sat) then
      call mem_allocate(this%nodesontop, this%maxbound, 'NODESONTOP', &
                        this%memoryPath)
    end if

    ! -- set input context pointers
    call mem_setptr(this%smassrate, 'SMASSRATE', this%input_mempath)
    !
    ! -- checkin input context pointers
    call mem_checkin(this%smassrate, 'SMASSRATE', this%memoryPath, &
                     'SMASSRATE', this%input_mempath)
    !
    ! -- Set values
    if (this%highest_sat) then
      do n = 1, this%maxbound
        this%nodesontop(n) = 0
      end do
    end if

  end subroutine src_allocate_arrays

  subroutine src_rp(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy
    class(GwtSrcType), intent(inout) :: this
    if (this%iper /= kper) return
    call this%BndExtType%bnd_rp()
    if (this%highest_sat) call this%set_nodesontop()
  end subroutine src_rp

  !> @brief Store nodelist in nodesontop
  !<
  subroutine set_nodesontop(this)
    implicit none
    ! -- dummy
    class(GwtSrcType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    ! !
    ! ! -- allocate if necessary
    ! if (.not. associated(this%nodesontop)) then
    !   allocate (this%nodesontop(this%maxbound))
    ! end if
    !
    ! -- copy nodelist into nodesontop
    do n = 1, this%nbound
      this%nodesontop(n) = this%nodelist(n)
    end do
  end subroutine set_nodesontop

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! This subroutine:
  !!   - calculates hcof and rhs terms
  !!   - skip if no sources
  !<
  subroutine src_cf(this)
    ! -- dummy
    class(GwtSrcType) :: this
    ! -- local
    integer(I4B) :: i, node
    real(DP) :: q
    !
    ! -- Return if no sources
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each source entry
    do i = 1, this%nbound
      !
      ! -- Find the node number
      if (this%highest_sat) then
        node = this%nodesontop(i)
      else
        node = this%nodelist(i)
      end if
      !
      ! -- reset nodelist to highest active
      if (this%highest_sat) then
        if (this%fmi%gwfsat(node) == 0) &
          call this%dis%highest_saturated(node, this%fmi%gwfsat)
        this%nodelist(i) = node
      end if

      this%hcof(i) = DZERO
      if (this%ibound(node) <= 0) then
        this%rhs(i) = DZERO
        cycle
      end if
      q = this%bound_value(1, i)
      this%rhs(i) = -q
    end do
  end subroutine src_cf

  !> @brief Add matrix terms related to specified mass source loading
  !!
  !! Copy rhs and hcof into solution rhs and amat
  !<
  subroutine src_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GwtSrcType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: i, n, ipos
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
      ! -- If mover is active and mass is being withdrawn,
      !    store available mass (as positive value).
      if (this%imover == 1 .and. this%rhs(i) > DZERO) then
        call this%pakmvrobj%accumulate_qformvr(i, this%rhs(i))
      end if
    end do
  end subroutine src_fc

  !> @brief Define list labels
  !!
  !! Define the list heading that is written to iout when PRINT_INPUT
  !! option is used.
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(GwtSrcType), intent(inout) :: this
    ! -- local
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
  !> @brief Support function for specified mass source loading observations
  !!
  !! This function:
  !!   - returns true because SRC package supports observations.
  !!   - overrides BndType%bnd_obs_supported()
  !<
  logical function src_obs_supported(this)
    implicit none
    ! -- dummy
    class(GwtSrcType) :: this
    !
    src_obs_supported = .true.
  end function src_obs_supported

  !> @brief Define observations
  !!
  !! This subroutine:
  !!   - stores observation types supported by SRC package.
  !!   - overrides BndType%bnd_df_obs
  !<
  subroutine src_df_obs(this)
    implicit none
    ! -- dummy
    class(GwtSrcType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType('src', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine src_df_obs

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index
  !!  and row.
  !<
  function src_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(GwtSrcType), intent(inout) :: this
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%smassrate(row)
    case default
    end select
  end function src_bound_value

end module GwtSrcModule
