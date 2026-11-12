!> @brief Package view of input parameter loading information
!!
!! Gives a package the loading information of its input fields (see
!! IdmParamMetaModule): the write stamps, to apply only the values its
!! input wrote in the current time step, the time series links of a
!! static block, and the read states of the period block load.
!<
module IdmParamMetaViewModule

  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENVARNAME
  use CharacterStringModule, only: CharacterStringType
  use MemoryManagerModule, only: mem_setptr, get_isize
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use LoadContextModule, only: ReadStateVarType, rsv_name
  use IdmParamMetaModule, only: PARAM_ISTEP, PARAM_VARS, InputStepType, &
                                tslink_id

  implicit none
  private
  public :: IdmParamMetaViewType

  !> @brief Loading information of a package's input fields
  !<
  type :: IdmParamMetaViewType
    character(len=LENVARNAME), dimension(:), allocatable :: fieldnames !< input field names
    integer(I4B), dimension(:), allocatable :: col !< stamp column of each field (0 if none)
    integer(I4B), dimension(:, :), pointer, contiguous :: istep => null() !< input step of each value's last write
    type(InputStepType) :: step !< current input step
    type(ReadStateVarType), dimension(:), allocatable :: reads !< read state of each field
    integer(I4B), dimension(:), allocatable :: lcol !< time series link column of each field (0 if none)
    integer(I4B), dimension(:, :), pointer, contiguous :: tslink => null() !< time series link flags
  contains
    procedure :: init
    procedure :: update
    procedure :: is_written
    procedure :: is_linked
    procedure :: read_state
    procedure :: destroy
    procedure, private :: tslink_init
  end type IdmParamMetaViewType

contains

  !> @brief Set pointers to the loading information of the named fields
  !!
  !! With stamps (default true), stores an error for a field without
  !! write stamps. With tslink_block, stores an error for a field without
  !! time series links in that block.
  !<
  subroutine init(this, fieldnames, input_mempath, stamps, tslink_block)
    class(IdmParamMetaViewType), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: fieldnames !< input field names
    character(len=*), intent(in) :: input_mempath !< input context memory path
    logical(LGP), optional, intent(in) :: stamps !< write stamps required
    character(len=*), optional, intent(in) :: tslink_block !< static block with time series links
    type(CharacterStringType), dimension(:), pointer, contiguous :: vars
    character(len=LENVARNAME) :: rsname
    logical(LGP) :: lstamps
    integer(I4B) :: j, m, nvars, isize

    lstamps = .true.
    if (present(stamps)) lstamps = stamps
    allocate (this%fieldnames(size(fieldnames)))
    allocate (this%col(size(fieldnames)))
    allocate (this%reads(size(fieldnames)))
    allocate (this%lcol(size(fieldnames)))
    this%col = 0
    this%lcol = 0
    do j = 1, size(fieldnames)
      nullify (this%reads(j)%invar)
    end do

    ! write stamps
    call get_isize(PARAM_VARS, input_mempath, nvars)
    if (nvars > 0) then
      call mem_setptr(this%istep, PARAM_ISTEP, input_mempath)
      call mem_setptr(vars, PARAM_VARS, input_mempath)
      call this%step%init()
    end if

    do j = 1, size(fieldnames)
      this%fieldnames(j) = fieldnames(j)
      if (nvars > 0) then
        do m = 1, size(vars)
          if (vars(m) == fieldnames(j)) this%col(j) = m
        end do
      end if
      if (lstamps .and. this%col(j) == 0) then
        write (errmsg, '(a)') 'Input write stamps not found for '// &
          trim(fieldnames(j))//'.'
        call store_error(errmsg)
      end if

      ! read state, if the loader keeps one
      rsname = rsv_name(fieldnames(j))
      call get_isize(rsname, input_mempath, isize)
      if (isize > 0) call mem_setptr(this%reads(j)%invar, rsname, input_mempath)
    end do

    ! time series links of a static block
    if (present(tslink_block)) then
      call this%tslink_init(fieldnames, input_mempath, tslink_block)
    end if
  end subroutine init

  !> @brief Set pointers to the time series links of a static block
  !<
  subroutine tslink_init(this, fieldnames, input_mempath, blockname)
    class(IdmParamMetaViewType), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: fieldnames !< input field names
    character(len=*), intent(in) :: input_mempath !< input context memory path
    character(len=*), intent(in) :: blockname !< static block name
    type(CharacterStringType), dimension(:), pointer, contiguous :: vars
    character(len=2) :: id
    integer(I4B) :: j, m, isize

    id = tslink_id(blockname)
    isize = 0
    if (id /= '') call get_isize('TSLINKVARS_'//id, input_mempath, isize)
    if (isize > 0) then
      call mem_setptr(this%tslink, 'TSLINK_'//id, input_mempath)
      call mem_setptr(vars, 'TSLINKVARS_'//id, input_mempath)
      do j = 1, size(fieldnames)
        do m = 1, size(vars)
          if (vars(m) == fieldnames(j)) this%lcol(j) = m
        end do
      end do
    end if
    do j = 1, size(fieldnames)
      if (this%lcol(j) == 0) then
        write (errmsg, '(a)') 'Time series links not found for '// &
          trim(fieldnames(j))//'.'
        call store_error(errmsg)
      end if
    end do
  end subroutine tslink_init

  !> @brief Bring the input step up to the current time step
  !!
  !! Call once in a time step before is_written.
  !<
  subroutine update(this)
    class(IdmParamMetaViewType), intent(inout) :: this
    if (associated(this%step%istep)) call this%step%update()
  end subroutine update

  !> @brief Return true if value n of field j was written in this time step
  !<
  function is_written(this, j, n) result(res)
    class(IdmParamMetaViewType), intent(in) :: this
    integer(I4B), intent(in) :: j !< field
    integer(I4B), intent(in) :: n !< value
    logical(LGP) :: res
    res = .false.
    if (this%col(j) > 0) res = this%istep(this%col(j), n) == this%step%istep
  end function is_written

  !> @brief Return true if value n of field j is linked to a time series
  !<
  function is_linked(this, j, n) result(res)
    class(IdmParamMetaViewType), intent(in) :: this
    integer(I4B), intent(in) :: j !< field
    integer(I4B), intent(in) :: n !< value (block row)
    logical(LGP) :: res
    res = .false.
    if (this%lcol(j) > 0) res = this%tslink(this%lcol(j), n) /= 0
  end function is_linked

  !> @brief Return the read state of field j, or -1 if it has none
  !<
  function read_state(this, j) result(rs)
    class(IdmParamMetaViewType), intent(in) :: this
    integer(I4B), intent(in) :: j !< field
    integer(I4B) :: rs
    rs = -1
    if (associated(this%reads(j)%invar)) rs = this%reads(j)%invar
  end function read_state

  !> @brief Release pointers and arrays
  !<
  subroutine destroy(this)
    class(IdmParamMetaViewType), intent(inout) :: this
    integer(I4B) :: j
    if (allocated(this%reads)) then
      do j = 1, size(this%reads)
        nullify (this%reads(j)%invar)
      end do
      deallocate (this%reads)
    end if
    if (allocated(this%fieldnames)) deallocate (this%fieldnames)
    if (allocated(this%col)) deallocate (this%col)
    if (allocated(this%lcol)) deallocate (this%lcol)
    nullify (this%istep)
    nullify (this%tslink)
    call this%step%destroy()
  end subroutine destroy

end module IdmParamMetaViewModule
