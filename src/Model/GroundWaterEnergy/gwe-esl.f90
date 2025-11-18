module GweEslModule
  !
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DEM1, DONE, LENFTYPE
  use SimVariablesModule, only: warnmsg
  use SimModule, only: store_warning
  use BndExtModule, only: BndExtType
  use ObsModule, only: DefaultObsIdProcessor
  use GweInputDataModule, only: GweInputDataType
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: esl_create
  !
  character(len=LENFTYPE) :: ftype = 'ESL'
  character(len=16) :: text = '             ESL'
  !
  type, extends(BndExtType) :: GweEslType

    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in mst
    real(DP), dimension(:), pointer, contiguous :: senerrate => null() !< energy source loading rate

  contains

    procedure :: allocate_scalars => esl_allocate_scalars
    procedure :: allocate_arrays => esl_allocate_arrays
    procedure :: bnd_cf => esl_cf
    procedure :: bnd_ck => esl_ck
    procedure :: bnd_fc => esl_fc
    procedure :: bnd_da => esl_da
    procedure :: define_listlabel
    procedure :: bound_value => esl_bound_value
    procedure :: ener_mult
    ! -- methods for observations
    procedure, public :: bnd_obs_supported => esl_obs_supported
    procedure, public :: bnd_df_obs => esl_df_obs

  end type GweEslType

contains

  !> @brief Create an energy source loading package
  !!
  !! This subroutine points bndobj to the newly created package
  !<
  subroutine esl_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        gwecommon, input_mempath)
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
    type(GweInputDataType), intent(in), target :: gwecommon !< shared data container for use by multiple GWE packages
    character(len=*), intent(in) :: input_mempath
    ! -- local
    type(GweEslType), pointer :: eslobj
    !
    ! -- Allocate the object and assign values to object variables
    allocate (eslobj)
    packobj => eslobj
    !
    ! -- Create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, input_mempath)
    packobj%text = text
    !
    ! -- Allocate scalars
    call eslobj%allocate_scalars()
    !
    ! -- Initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1
    !
    ! -- Store pointer to shared data module for accessing cpw, rhow
    !    for the budget calculations, and for accessing the latent heat of
    !    vaporization for evaporative cooling.
    eslobj%gwecommon => gwecommon
  end subroutine esl_create

  !> @brief Deallocate memory
  !<
  subroutine esl_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GweEslType) :: this
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
  end subroutine esl_da

  !> @brief Allocate scalars
  !!
  !! Allocate scalars specific to this energy source loading package
  !<
  subroutine esl_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweEslType) :: this
    !
    ! -- base class allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    !
    ! -- Set values
  end subroutine esl_allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine esl_allocate_arrays(this, nodelist, auxvar)
    use MemoryManagerModule, only: mem_setptr, mem_checkin
    ! -- dummy
    class(GweEslType) :: this
    ! -- local
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist !< package nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar !< package aux variable array

    ! -- base class allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)

    ! -- set input context pointers
    call mem_setptr(this%senerrate, 'SENERRATE', this%input_mempath)
    !
    ! -- checkin input context pointers
    call mem_checkin(this%senerrate, 'SENERRATE', this%memoryPath, &
                     'SENERRATE', this%input_mempath)
  end subroutine esl_allocate_arrays

  !> @brief Check energy source loading boundary condition data
  !<
  subroutine esl_ck(this)
    ! -- dummy
    class(GweEslType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: node
    ! -- formats
    character(len=*), parameter :: fmtenermulterr = &
      "('ESL BOUNDARY (',i0,') ESL MULTIPLIER (',g10.3,') IS &
      &LESS THAN ZERO THEREBY REVERSING THE ORIGINAL SIGN ON THE &
      &AMOUNT OF ENERGY ENTERING OR EXITING THE MODEL.')"
    !
    ! -- check stress period data
    do i = 1, this%nbound
      node = this%nodelist(i)
      !
      ! -- accumulate warnings
      if (this%iauxmultcol > 0) then
        if (this%auxvar(this%iauxmultcol, i) < DZERO) then
          write (warnmsg, fmt=fmtenermulterr) &
            i, this%auxvar(this%iauxmultcol, i)
          call store_warning(warnmsg)
          write (this%iout, '(/1x,a)') 'WARNING: '//trim(warnmsg)
        end if
      end if
    end do
  end subroutine esl_ck

  !> @brief Formulate the HCOF and RHS terms
  !!
  !! This subroutine:
  !!   - calculates hcof and rhs terms
  !!   - skip if no sources
  !<
  subroutine esl_cf(this)
    ! -- dummy
    class(GweEslType) :: this
    ! -- local
    integer(I4B) :: i, node
    real(DP) :: q
    !
    ! -- Return if no sources
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each source entry
    do i = 1, this%nbound
      node = this%nodelist(i)
      this%hcof(i) = DZERO
      if (this%ibound(node) <= 0) then
        this%rhs(i) = DZERO
        cycle
      end if
      !
      ! -- set energy loading rate accounting for multiplier
      q = this%ener_mult(i)
      !
      this%rhs(i) = -q
    end do
  end subroutine esl_cf

  !> @brief Add matrix terms related to specified energy source loading
  !!
  !! Copy rhs and hcof into solution rhs and amat
  !<
  subroutine esl_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GweEslType) :: this
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
  end subroutine esl_fc

  !> @brief Define list labels
  !!
  !! Define the list heading that is written to iout when
  !! PRINT_INPUT option is used.
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(GweEslType), intent(inout) :: this
    !
    ! -- Create the header list label
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

  !> @brief Support function for specified energy source loading observations
  !!
  !! This function:
  !!   - returns true because ESL package supports observations.
  !!   - overrides BndType%bnd_obs_supported()
  !<
  logical function esl_obs_supported(this)
    implicit none
    ! -- dummy
    class(GweEslType) :: this
    !
    esl_obs_supported = .true.
  end function esl_obs_supported

  !> @brief Define observations
  !!
  !! This subroutine:
  !!   - stores observation types supported by ESL package.
  !!   - overrides BndType%bnd_df_obs
  !<
  subroutine esl_df_obs(this)
    implicit none
    ! -- dummy
    class(GweEslType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType('esl', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine esl_df_obs

  !> @ brief Return a bound value
  !!
  !!  Return a bound value associated with an ncolbnd index
  !!  and row.
  !<
  function esl_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(GweEslType), intent(inout) :: this
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    !
    select case (col)
    case (1)
      bndval = this%senerrate(row)
    case default
    end select
  end function esl_bound_value

  !> @brief Return a value that applies a multiplier
  !!
  !! Apply multiplier to specified energy load depending on user-selected
  !! option
  !<
  function ener_mult(this, row) result(ener)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(GweEslType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: ener
    !
    if (this%iauxmultcol > 0) then
      ener = this%senerrate(row) * this%auxvar(this%iauxmultcol, row)
    else
      ener = this%senerrate(row)
    end if
  end function ener_mult

end module GweEslModule
