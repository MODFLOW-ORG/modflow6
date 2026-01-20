!> @brief This module contains the PackageBudgetModule Module
!!
!! The PackageBudgetType object defined here provides flows to the GWT
!! model.  The PackageBudgetType can be filled with flows from a budget
!! object that was written from a previous GWF simulation, or its
!! individual members can be pointed to flows that are being calculated
!! by a GWF model that is running as part of this simulation.
!<
module PackageBudgetModule

  use KindModule
  use ConstantsModule, only: LENPACKAGENAME, LENAUXNAME, LENMEMPATH
  use MemoryManagerModule, only: mem_allocate, mem_reassignptr, &
                                 mem_reallocate, mem_deallocate

  implicit none

  private
  public :: PackageBudgetType

  !> @brief Derived type for storing flows
  !!
  !! This derived type stores flows and provides them through the FMI
  !! package to other parts of GWT.
  !!
  !<
  type :: PackageBudgetType

    character(len=LENMEMPATH) :: memoryPath = '' !< the location in the memory manager where the variables are stored
    character(len=LENPACKAGENAME), pointer :: name => null() !< name of the package
    character(len=LENPACKAGENAME), pointer :: budtxt => null() !< type of flow (CHD, RCH, RCHA, ...)
    character(len=LENAUXNAME), dimension(:), pointer, &
      contiguous :: auxname => null() !< vector of auxname
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    integer(I4B), pointer :: nbound => null() !< number of boundaries for current stress period
    integer(I4B), dimension(:), pointer, contiguous :: nodelist => null() !< vector of reduced node numbers
    real(DP), dimension(:), pointer, contiguous :: flow => null() !< calculated flow
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< auxiliary variable array

  contains

    procedure :: initialize
    procedure :: set_name
    procedure :: set_auxname
    procedure :: set_pointers
    procedure :: copy_values
    procedure :: get_flow
    procedure :: da

  end type PackageBudgetType

contains

  !> @ brief Initialize a PackageBudgetType object
  !!
  !!  Establish the memory path and allocate and initialize member variables.
  !!
  !<
  subroutine initialize(this, mempath)
    class(PackageBudgetType) :: this !< PackageBudgetType object
    character(len=*), intent(in) :: mempath !< memory path in memory manager
    this%memoryPath = mempath
    !
    ! -- allocate member variables in memory manager
    call mem_allocate(this%name, LENPACKAGENAME, 'NAME', mempath)
    call mem_allocate(this%budtxt, LENPACKAGENAME, 'BUDTXT', mempath)
    call mem_allocate(this%naux, 'NAUX', mempath)
    call mem_allocate(this%auxname, LENAUXNAME, 0, 'AUXNAME', this%memoryPath)
    call mem_allocate(this%nbound, 'NBOUND', mempath)
    call mem_allocate(this%nodelist, 0, 'NODELIST', mempath)
    call mem_allocate(this%flow, 0, 'FLOW', mempath)
    call mem_allocate(this%auxvar, 0, 0, 'AUXVAR', mempath)
    !
    ! -- initialize
    this%name = ''
    this%budtxt = ''
    this%naux = 0
    this%nbound = 0
  end subroutine initialize

  !> @ brief Set names for this PackageBudgetType object
  !!
  !!  Set the name of the package and the name of the of budget text
  !!
  !<
  subroutine set_name(this, name, budtxt)
    class(PackageBudgetType) :: this !< PackageBudgetType object
    character(len=LENPACKAGENAME) :: name !< name of the package (WEL-1, DRN-4, etc.)
    character(len=LENPACKAGENAME) :: budtxt !< name of budget term (CHD, RCH, EVT, DRN-TO-MVR, etc.)
    this%name = name
    this%budtxt = budtxt
  end subroutine set_name

  !> @ brief Set aux names for this PackageBudgetType object
  !!
  !!  Set the number of auxiliary variables and the names of the
  !!  auxiliary variables
  !!
  !<
  subroutine set_auxname(this, naux, auxname)
    class(PackageBudgetType) :: this !< PackageBudgetType object
    integer(I4B), intent(in) :: naux !< number of auxiliary variables
    character(len=LENAUXNAME), contiguous, &
      dimension(:), intent(in) :: auxname !< array of names for auxiliary variables
    this%naux = naux
    call mem_reallocate(this%auxname, LENAUXNAME, naux, 'AUXNAME', &
                        this%memoryPath)
    this%auxname(:) = auxname(:)
  end subroutine set_auxname

  !> @ brief Point members of this class to data stored in GWF packages
  !!
  !!  The routine is called when a GWF model is being run concurrently with
  !!  a GWT model.  In this situation, the member variables NBOUND, NODELIST,
  !!  FLOW, and AUXVAR are pointed into member variables of the individual
  !!  GWF Package members stored in BndType.
  !!
  !<
  subroutine set_pointers(this, flowvarname, mem_path_target, input_mempath)
    use ConstantsModule, only: LENVARNAME
    class(PackageBudgetType) :: this !< PackageBudgetType object
    character(len=*), intent(in) :: flowvarname !< name of variable storing flow (SIMVALS, SIMTOMVR)
    character(len=*), intent(in) :: mem_path_target !< path where target variable is stored
    character(len=*), intent(in) :: input_mempath
    character(len=LENVARNAME) :: auxvarname
    !
    ! -- set memory manager aux varname
    if (input_mempath /= '') then
      auxvarname = 'AUXVAR_IDM'
    else
      auxvarname = 'AUXVAR'
    end if
    !
    ! -- Reassign pointers to variables in the flow model
    call mem_reassignptr(this%nbound, 'NBOUND', this%memoryPath, &
                         'NBOUND', mem_path_target)
    call mem_reassignptr(this%nodelist, 'NODELIST', this%memoryPath, &
                         'NODELIST', mem_path_target)
    call mem_reassignptr(this%flow, 'FLOW', this%memoryPath, &
                         flowvarname, mem_path_target)
    call mem_reassignptr(this%auxvar, 'AUXVAR', this%memoryPath, &
                         auxvarname, mem_path_target)
  end subroutine set_pointers

  !> @ brief Copy data read from a budget file into this object
  !!
  !!  The routine is called when GWF flows are read from a budget file
  !!  created from a previous GWF simulation.  Arrays here must be
  !!  dynamically resized because maxbound is not known unless the
  !!  entire budget file was read first.
  !!
  !<
  subroutine copy_values(this, nbound, nodelist, flow, auxvar)
    class(PackageBudgetType) :: this !< PackageBudgetType object
    integer(I4B), intent(in) :: nbound !< number of entries
    integer(I4B), dimension(:), contiguous, intent(in) :: nodelist !< array of GWT node numbers
    real(DP), dimension(:), contiguous, intent(in) :: flow !< array of flow rates
    real(DP), dimension(:, :), contiguous, intent(in) :: auxvar !< array of auxiliary variables
    integer(I4B) :: i
    !
    ! -- Assign variables
    this%nbound = nbound
    !
    ! -- Lists are not large enough (maxbound is not known), so need to
    !    reallocate based on size in binary budget file.
    if (size(this%nodelist) < nbound) then
      call mem_reallocate(this%nodelist, nbound, 'NODELIST', this%memoryPath)
      call mem_reallocate(this%flow, nbound, 'FLOW', this%memoryPath)
      call mem_reallocate(this%auxvar, this%naux, nbound, 'AUXVAR', &
                          this%memoryPath)
    end if
    !
    ! -- Copy values into member variables
    do i = 1, nbound
      this%nodelist(i) = nodelist(i)
      this%flow(i) = flow(i)
      this%auxvar(:, i) = auxvar(:, i)
    end do
  end subroutine copy_values

  !> @ brief Get flow rate for specified entry
  !!
  !!  Return the flow rate for the specified entry
  !!
  !<
  function get_flow(this, i) result(flow)
    class(PackageBudgetType) :: this !< PackageBudgetType object
    integer(I4B), intent(in) :: i !< entry number
    real(DP) :: flow
    flow = this%flow(i)
  end function get_flow

  !> @ brief Deallocate
  !!
  !!  Free any memory associated with this object
  !!
  !<
  subroutine da(this)
    class(PackageBudgetType) :: this
    call mem_deallocate(this%name, 'NAME', this%memoryPath)
    call mem_deallocate(this%budtxt, 'BUDTXT', this%memoryPath)
    call mem_deallocate(this%naux)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%nbound)
    call mem_deallocate(this%nodelist, 'NODELIST', this%memoryPath)
    call mem_deallocate(this%flow, 'FLOW', this%memoryPath)
    call mem_deallocate(this%auxvar, 'AUXVAR', this%memoryPath)
  end subroutine da

end module PackageBudgetModule
