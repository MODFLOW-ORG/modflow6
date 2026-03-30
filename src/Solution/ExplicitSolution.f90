!> @brief Explicit solutions for solving explicit models.
!!
!! Explicit solutions manage and solve explicit models. An
!! explicit model solves itself, as opposed to a numerical
!! model which requires a numerical solution procedure. An
!! explicit solution involves a double loop: an outer loop
!! that continues until no model has pending work, as well
!! as an inner loop that scrolls through the models in the
!! solution and tells each to solve itself. The outer loop
!! is necessary because the explicit models may be coupled
!! by exchanges; models may send work during each solve to
!! other models. The outer loop continues until all models
!! have completed all work.
!<
module ExplicitSolutionModule
  use KindModule, only: I4B, DP
  use TimerModule, only: code_timer
  use ConstantsModule, only: LENMEMPATH, LENSOLUTIONNAME, MVALIDATE, &
                             MNORMAL, LINELENGTH, DZERO
  use MemoryHelperModule, only: create_mem_path
  use BaseModelModule, only: BaseModelType
  use ExplicitModelModule, only: ExplicitModelType, &
                                 AddExplicitModelToList, &
                                 GetExplicitModelFromList
  use BaseExchangeModule, only: BaseExchangeType
  use BaseSolutionModule, only: BaseSolutionType, AddBaseSolutionToList
  use ListModule, only: ListType
  use ListsModule, only: basesolutionlist
  use SimVariablesModule, only: iout, isim_mode
  use BlockParserModule, only: BlockParserType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use InputOutputModule, only: getunit

  implicit none
  private

  public :: create_explicit_solution
  public :: ExplicitSolutionType

  !> @brief Manages and solves explicit models.
  type, extends(BaseSolutionType) :: ExplicitSolutionType
    character(len=LENMEMPATH) :: memoryPath !< the path for storing solution variables in the memory manager
    type(ListType), pointer :: modellist !< list of models in solution
    type(ListType), pointer :: exchangelist !< list of exchanges in solution
    integer(I4B), pointer :: id !< solution number
    integer(I4B), pointer :: iu !< input file unit
    real(DP), pointer :: ttsoln !< timer - total solution time
    integer(I4B), pointer :: icnvg => null() !< convergence flag
    type(BlockParserType) :: parser !< block parser object
  contains
    procedure :: sln_df
    procedure :: sln_ar
    procedure :: sln_dt
    procedure :: sln_ad
    procedure :: sln_ot
    procedure :: sln_ca
    procedure :: sln_fp
    procedure :: sln_da
    procedure :: add_model
    procedure :: add_exchange
    procedure :: get_models
    procedure :: get_exchanges
    procedure :: save

    procedure, private :: allocate_scalars

    ! Expose these for use through the BMI/XMI:
    procedure, public :: prepareSolve
    procedure, public :: solve
    procedure, public :: finalizeSolve

  end type ExplicitSolutionType

contains

  !> @ brief Create a new solution
  !!
  !! Create a new solution using the data in filename, assign this new
  !! solution an id number and store the solution in the basesolutionlist.
  !! Also open the filename for later reading.
  !<
  subroutine create_explicit_solution(exp_sol, filename, id)
    ! modules
    use InputOutputModule, only: getunit, openfile
    ! dummy
    class(ExplicitSolutionType), pointer :: exp_sol !< the create solution
    character(len=*), intent(in) :: filename !< solution input file name
    integer(I4B), intent(in) :: id !< solution id
    ! local
    integer(I4B) :: inunit
    class(BaseSolutionType), pointer :: solbase => null()
    character(len=LENSOLUTIONNAME) :: solutionname

    ! Create a new solution and add it to the basesolutionlist container
    solbase => exp_sol
    write (solutionname, '(a, i0)') 'SLN_', id
    exp_sol%name = solutionname
    exp_sol%memoryPath = create_mem_path(solutionname)
    allocate (exp_sol%modellist)
    allocate (exp_sol%exchangelist)
    call exp_sol%allocate_scalars()
    call AddBaseSolutionToList(basesolutionlist, solbase)
    exp_sol%id = id

    ! Open solution input file for reading later after problem size is known
    ! Check to see if the file is already opened, which can happen when
    ! running in single model mode
    inquire (file=filename, number=inunit)
    if (inunit < 0) inunit = getunit()
    exp_sol%iu = inunit
    write (iout, '(/a,a/)') ' Creating explicit solution (EMS): ', exp_sol%name
    call openfile(exp_sol%iu, iout, filename, 'EMS')

    ! Initialize block parser
    call exp_sol%parser%Initialize(exp_sol%iu, iout)
  end subroutine create_explicit_solution

  !> @ brief Allocate scalars
  subroutine allocate_scalars(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance

    ! allocate
    call mem_allocate(this%id, 'ID', this%memoryPath)
    call mem_allocate(this%iu, 'IU', this%memoryPath)
    call mem_allocate(this%ttsoln, 'TTSOLN', this%memoryPath)
    call mem_allocate(this%icnvg, 'ICNVG', this%memoryPath)

    ! initialize
    this%id = 0
    this%iu = 0
    this%ttsoln = DZERO
    this%icnvg = 0
  end subroutine allocate_scalars

  !> @ brief Define the solution
  subroutine sln_df(this)
    class(ExplicitSolutionType) :: this
  end subroutine

  !> @ brief Allocate and read
  subroutine sln_ar(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    ! close ems input file
    call this%parser%Clear()
  end subroutine sln_ar

  !> @ brief Calculate time step length
  subroutine sln_dt(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
  end subroutine sln_dt

  !> @ brief Advance the solution
  subroutine sln_ad(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance

    ! reset convergence flag
    this%icnvg = 0
  end subroutine sln_ad

  !> @ brief Output
  subroutine sln_ot(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
  end subroutine sln_ot

  subroutine sln_fp(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
  end subroutine sln_fp

  !> @ brief Deallocate
  subroutine sln_da(this)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance

    ! lists
    call this%modellist%Clear()
    deallocate (this%modellist)
    call this%exchangelist%Clear()
    deallocate (this%exchangelist)

    ! scalars
    call mem_deallocate(this%id)
    call mem_deallocate(this%iu)
    call mem_deallocate(this%ttsoln)
    call mem_deallocate(this%icnvg)
  end subroutine sln_da

  !> @ brief Calculate
  subroutine sln_ca(this, isgcnvg, isuppress_output)
    ! dummy
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    integer(I4B), intent(inout) :: isgcnvg !< solution group convergence flag
    integer(I4B), intent(in) :: isuppress_output !< flag for suppressing output
    ! local
    class(ExplicitModelType), pointer :: mp => null()
    character(len=LINELENGTH) :: line
    character(len=LINELENGTH) :: fmt
    integer(I4B) :: im
    integer(I4B) :: kiter

    kiter = 1

    call this%prepareSolve()

    select case (isim_mode)
    case (MVALIDATE)
      line = 'mode="validation" -- Skipping assembly and solution.'
      fmt = "(/,1x,a,/)"
      do im = 1, this%modellist%Count()
        mp => GetExplicitModelFromList(this%modellist, im)
        call mp%model_message(line, fmt=fmt)
      end do
    case (MNORMAL)

      ! solve the models
      call this%solve(kiter)

      ! finish up
      call this%finalizeSolve(kiter, isgcnvg, isuppress_output)
    end select
  end subroutine sln_ca

  !> @ brief Prepare to solve
  subroutine prepareSolve(this)
    ! dummy
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    ! local
    integer(I4B) :: im
    class(ExplicitModelType), pointer :: mp => null()

    ! advance model
    do im = 1, this%modellist%Count()
      mp => GetExplicitModelFromList(this%modellist, im)
      call mp%model_ad()
    end do

    ! advance solution
    call this%sln_ad()
  end subroutine prepareSolve

  !> @ brief Solve models
  subroutine solve(this, kiter)
    ! dummy
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    integer(I4B), intent(in) :: kiter !< Picard iteration (1 for explicit)
    ! local
    class(ExplicitModelType), pointer :: mp => null()
    integer(I4B) :: im
    logical :: any_pending
    real(DP) :: ttsoln

    call code_timer(0, ttsoln, this%ttsoln)

    ! Outer loop: repeat until no model has pending work
    do
      ! Solve every model in the solution
      do im = 1, this%modellist%Count()
        mp => GetExplicitModelFromList(this%modellist, im)
        call mp%model_solve()
      end do

      ! Continue if any have pending work
      any_pending = .false.
      do im = 1, this%modellist%Count()
        mp => GetExplicitModelFromList(this%modellist, im)
        if (mp%has_pending()) then
          any_pending = .true.
          exit
        end if
      end do
      if (.not. any_pending) exit
    end do

    call code_timer(1, ttsoln, this%ttsoln)
    this%icnvg = 1
  end subroutine solve

  !> @ brief Finalize solve
  subroutine finalizeSolve(this, kiter, isgcnvg, isuppress_output)
    ! dummy
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    integer(I4B), intent(in) :: kiter !< Picard iteration number (always 1 for explicit)
    integer(I4B), intent(inout) :: isgcnvg !< solution group convergence flag
    integer(I4B), intent(in) :: isuppress_output !< flag for suppressing output
    ! local
    integer(I4B) :: im
    class(ExplicitModelType), pointer :: mp => null()

    ! Calculate flow for each model
    do im = 1, this%modellist%Count()
      mp => GetExplicitModelFromList(this%modellist, im)
      call mp%model_cq(this%icnvg, isuppress_output)
    end do

    ! Budget terms for each model
    do im = 1, this%modellist%Count()
      mp => GetExplicitModelFromList(this%modellist, im)
      call mp%model_bd(this%icnvg, isuppress_output)
    end do
  end subroutine finalizeSolve

  !> @ brief Save output
  subroutine save(this, filename)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    character(len=*), intent(in) :: filename !< filename to save solution data
    integer(I4B) :: inunit
    inunit = getunit()
    open (unit=inunit, file=filename, status='unknown')
    write (inunit, *) 'The save routine currently writes nothing'
    close (inunit)
  end subroutine save

  !> @ brief Add explicit model to list
  subroutine add_model(this, mp)
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    class(BaseModelType), pointer, intent(in) :: mp !< model instance
    class(ExplicitModelType), pointer :: m => null()
    select type (mp)
    class is (ExplicitModelType)
      m => mp
      call AddExplicitModelToList(this%modellist, m)
    end select
  end subroutine add_model

  !> @brief Get a pointer to a list of models in the solution
  function get_models(this) result(models)
    type(ListType), pointer :: models !< pointer to the model list
    class(ExplicitSolutionType) :: this !< ExplicitSolutionType instance
    models => this%modellist
  end function get_models

  !> @ brief Add exchange to list of exchanges
  subroutine add_exchange(this, exchange)
    class(ExplicitSolutionType) :: this
    class(BaseExchangeType), pointer, intent(in) :: exchange
    class(*), pointer :: obj
    obj => exchange
    call this%exchangelist%Add(obj)
  end subroutine add_exchange

  !> @ brief Get list of exchanges
  function get_exchanges(this) result(exchanges)
    class(ExplicitSolutionType) :: this
    type(ListType), pointer :: exchanges
    exchanges => this%exchangelist
  end function get_exchanges

end module ExplicitSolutionModule
