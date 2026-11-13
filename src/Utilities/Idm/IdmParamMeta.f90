!> @brief Input parameter loading information published for packages
!!
!! Loaders record how and when input parameter values are written into
!! the input context. Packages that must apply a value only when it is
!! written read this information (see IdmParamMetaViewModule):
!!  - TSLINK_\<id\>, TSLINKVARS_\<id\> (block loader; PACKAGEDATA and
!!    OUTLETS blocks of advanced packages): (ncol, nrow), nonzero where a
!!    row's value is linked to a time series, and the column names. A
!!    linked value is written each time step.
!!  - PARAM_ISTEP, PARAM_RS, PARAM_VARS (keystring loader; PERIOD block
!!    parameters addressed by node): (ncol, nfeatures), the input step of
!!    each value's last write and how it was written (RS_VALUES, given in
!!    a block; RS_SERIES, by a time series), and the column names.
!!  - IN\<param\> (layer and grid array loaders; PERIOD block array
!!    parameters): the read state (RS_*) of the most recent period block.
!!  - ISTEP (simulation input path): the input step, the number of
!!    distinct time steps (kper, kstp) reached so far. It advances when a
!!    new time step is first seen, so a retried time step keeps it.
!<
module IdmParamMetaModule

  use KindModule, only: I4B
  use ConstantsModule, only: LENVARNAME, LENMEMPATH
  use CharacterStringModule, only: CharacterStringType
  use MemoryManagerModule, only: mem_allocate, mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use SimVariablesModule, only: idm_context

  implicit none
  private
  public :: tslink_id
  public :: tslink_allocate
  public :: input_step_allocate
  public :: input_step_update
  public :: InputStepType
  public :: ParamStampsType
  public :: PARAM_ISTEP, PARAM_RS, PARAM_VARS
  public :: ReadStateVarType
  public :: RS_UNSET, RS_NONE, RS_VALUES, RS_SERIES, RS_NETCDF
  public :: rsv_name

  character(len=*), parameter :: INPUT_ISTEP = 'ISTEP' !< input step
  character(len=*), parameter :: INPUT_KPER = 'ISTEP_KPER' !< stress period of the input step
  character(len=*), parameter :: INPUT_KSTP = 'ISTEP_KSTP' !< time step of the input step
  character(len=*), parameter :: PARAM_ISTEP = 'PARAM_ISTEP' !< write stamps
  character(len=*), parameter :: PARAM_RS = 'PARAM_RS' !< write stamp read states
  character(len=*), parameter :: PARAM_VARS = 'PARAM_VARS' !< write stamp column names

  !> @brief Read state of a period array (e.g. INRECHARGE)
  !!
  !! Set for each period block load. Write stamps (PARAM_RS) use RS_VALUES
  !! for a value given in a block and RS_SERIES for a value written by a
  !! time series.
  !<
  enum, bind(C)
    enumerator :: RS_UNSET = -1 !< not in use
    enumerator :: RS_NONE = 0 !< not given in the most recent period block
    enumerator :: RS_VALUES = 1 !< values given in the most recent period block
    enumerator :: RS_SERIES = 2 !< time-array series given in the most recent period block
    enumerator :: RS_NETCDF = 3 !< values read from NetCDF in the most recent period block
  end enum

  !> @brief Pointer type for read state variable
  !<
  type ReadStateVarType
    integer(I4B), pointer :: invar
  end type ReadStateVarType

  !> @brief The input step, kept current with the TDIS time step
  !<
  type :: InputStepType
    integer(I4B), pointer :: istep => null() !< input step
    integer(I4B), pointer :: kper => null() !< stress period of the input step
    integer(I4B), pointer :: kstp => null() !< time step of the input step
  contains
    procedure :: init => step_init
    procedure :: update => step_update
    procedure :: destroy => step_destroy
  end type InputStepType

  !> @brief Write stamps of a loader's input parameter values
  !<
  type :: ParamStampsType
    integer(I4B), dimension(:, :), pointer, contiguous :: istep => null() !< input step of each value's last write
    integer(I4B), dimension(:, :), pointer, contiguous :: rs => null() !< read state of each value's last write
    type(InputStepType) :: step !< current input step
  contains
    procedure :: allocate => stamps_allocate
    procedure :: set => stamps_set
    procedure :: destroy => stamps_destroy
  end type ParamStampsType

contains

  !> @brief Read state variable name of a parameter (e.g. INRECHARGE)
  !<
  function rsv_name(mf6varname) result(varname)
    character(len=*), intent(in) :: mf6varname
    character(len=LENVARNAME) :: varname
    integer(I4B) :: ilen
    character(len=2) :: prefix = 'IN'
    ilen = len_trim(mf6varname)
    if (ilen > (LENVARNAME - len(prefix))) then
      varname = prefix//mf6varname(1:(LENVARNAME - len(prefix)))
    else
      varname = prefix//trim(mf6varname)
    end if
  end function rsv_name

  !> @brief Time series link identifier for a block
  !!
  !! Return the identifier used to name a block's time series link
  !! variables (TSLINK_\<id\>, TSLINKVARS_\<id\>), or an empty string if the
  !! block has none.
  !<
  function tslink_id(blockname) result(id)
    character(len=*), intent(in) :: blockname
    character(len=2) :: id
    select case (blockname)
    case ('PACKAGEDATA')
      id = 'PD'
    case ('OUTLETS')
      id = 'OL'
    case default
      id = ''
    end select
  end function tslink_id

  !> @brief Allocate time series link flags for a block
  !!
  !! Allocate TSLINK_\<id\> (zero) and TSLINKVARS_\<id\> (set to varnames) for
  !! ncol columns and nrow rows. The caller sets the linked rows.
  !<
  subroutine tslink_allocate(id, varnames, nrow, mempath, tslink)
    character(len=*), intent(in) :: id !< block identifier from tslink_id
    character(len=*), dimension(:), intent(in) :: varnames !< column names
    integer(I4B), intent(in) :: nrow !< number of rows
    character(len=*), intent(in) :: mempath !< input context memory path
    integer(I4B), dimension(:, :), pointer, contiguous, intent(out) :: tslink !< link flags
    type(CharacterStringType), dimension(:), pointer, contiguous :: tslinkvars
    integer(I4B) :: m

    call mem_allocate(tslink, size(varnames), nrow, 'TSLINK_'//trim(id), &
                      mempath)
    call mem_allocate(tslinkvars, LENVARNAME, size(varnames), &
                      'TSLINKVARS_'//trim(id), mempath)
    tslink = 0
    do m = 1, size(varnames)
      tslinkvars(m) = varnames(m)
    end do
  end subroutine tslink_allocate

  !> @brief Memory path of the input step
  !<
  function input_step_mempath() result(mempath)
    character(len=LENMEMPATH) :: mempath
    mempath = create_mem_path(component='SIM', context=idm_context)
  end function input_step_mempath

  !> @brief Allocate the input step
  !!
  !! The input step starts at zero, with no time step reached.
  !<
  subroutine input_step_allocate()
    integer(I4B), pointer :: ival
    call mem_allocate(ival, INPUT_ISTEP, input_step_mempath())
    ival = 0
    call mem_allocate(ival, INPUT_KPER, input_step_mempath())
    ival = -1
    call mem_allocate(ival, INPUT_KSTP, input_step_mempath())
    ival = -1
  end subroutine input_step_allocate

  !> @brief Bring the input step up to the current time step
  !<
  subroutine input_step_update()
    type(InputStepType) :: step
    call step%init()
    call step%update()
    call step%destroy()
  end subroutine input_step_update

  !> @brief Set pointers to the input step
  !<
  subroutine step_init(this)
    class(InputStepType), intent(inout) :: this
    call mem_setptr(this%istep, INPUT_ISTEP, input_step_mempath())
    call mem_setptr(this%kper, INPUT_KPER, input_step_mempath())
    call mem_setptr(this%kstp, INPUT_KSTP, input_step_mempath())
  end subroutine step_init

  !> @brief Advance the input step if the time step is new
  !<
  subroutine step_update(this)
    use TdisModule, only: kper, kstp
    class(InputStepType), intent(inout) :: this
    if (kper /= this%kper .or. kstp /= this%kstp) then
      this%istep = this%istep + 1
      this%kper = kper
      this%kstp = kstp
    end if
  end subroutine step_update

  !> @brief Release pointers
  !<
  subroutine step_destroy(this)
    class(InputStepType), intent(inout) :: this
    nullify (this%istep)
    nullify (this%kper)
    nullify (this%kstp)
  end subroutine step_destroy

  !> @brief Allocate write stamps for ncol columns and nfeatures values
  !<
  subroutine stamps_allocate(this, varnames, nfeatures, mempath)
    class(ParamStampsType), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: varnames !< column names
    integer(I4B), intent(in) :: nfeatures !< number of values per column
    character(len=*), intent(in) :: mempath !< input context memory path
    type(CharacterStringType), dimension(:), pointer, contiguous :: vars
    integer(I4B) :: m

    call this%step%init()
    call mem_allocate(this%istep, size(varnames), nfeatures, PARAM_ISTEP, &
                      mempath)
    call mem_allocate(this%rs, size(varnames), nfeatures, PARAM_RS, mempath)
    call mem_allocate(vars, LENVARNAME, size(varnames), PARAM_VARS, mempath)
    this%istep = 0
    this%rs = RS_UNSET
    do m = 1, size(varnames)
      vars(m) = varnames(m)
    end do
  end subroutine stamps_allocate

  !> @brief Stamp value n of column col as written in this input step
  !<
  subroutine stamps_set(this, col, n, rs)
    class(ParamStampsType), intent(inout) :: this
    integer(I4B), intent(in) :: col !< column
    integer(I4B), intent(in) :: n !< value
    integer(I4B), intent(in) :: rs !< how the value was written (RS_VALUES or RS_SERIES)
    call this%step%update()
    this%istep(col, n) = this%step%istep
    this%rs(col, n) = rs
  end subroutine stamps_set

  !> @brief Release pointers; the arrays are freed with the input context
  !<
  subroutine stamps_destroy(this)
    class(ParamStampsType), intent(inout) :: this
    nullify (this%istep)
    nullify (this%rs)
    call this%step%destroy()
  end subroutine stamps_destroy

end module IdmParamMetaModule
