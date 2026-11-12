!> @brief Package view of input parameter loading information
!!
!! Gives a package the loading information of its input fields (see
!! IdmParamMetaModule): the write stamps, to apply only the values its
!! input wrote in the current time step, and the read states of the
!! period block load.
!<
module IdmParamMetaViewModule

  use KindModule, only: I4B, LGP
  use ConstantsModule, only: LENVARNAME
  use CharacterStringModule, only: CharacterStringType
  use MemoryManagerModule, only: mem_setptr, get_isize
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use LoadContextModule, only: ReadStateVarType, rsv_name
  use IdmParamMetaModule, only: PARAM_ISTEP, PARAM_VARS, InputStepType

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
  contains
    procedure :: init
    procedure :: update
    procedure :: is_written
    procedure :: read_state
    procedure :: destroy
  end type IdmParamMetaViewType

contains

  !> @brief Set pointers to the loading information of the named fields
  !!
  !! With stamps (default true), stores an error for a field without
  !! write stamps.
  !<
  subroutine init(this, fieldnames, input_mempath, stamps)
    class(IdmParamMetaViewType), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: fieldnames !< input field names
    character(len=*), intent(in) :: input_mempath !< input context memory path
    logical(LGP), optional, intent(in) :: stamps !< write stamps required
    type(CharacterStringType), dimension(:), pointer, contiguous :: vars
    character(len=LENVARNAME) :: rsname
    logical(LGP) :: lstamps
    integer(I4B) :: j, m, nvars, isize

    lstamps = .true.
    if (present(stamps)) lstamps = stamps
    allocate (this%fieldnames(size(fieldnames)))
    allocate (this%col(size(fieldnames)))
    allocate (this%reads(size(fieldnames)))
    this%col = 0
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
  end subroutine init

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
    nullify (this%istep)
    call this%step%destroy()
  end subroutine destroy

end module IdmParamMetaViewModule
