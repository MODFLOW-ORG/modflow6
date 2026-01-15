module ProfilerModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DNODATA, DZERO, LENMEMPATH, LINELENGTH
  use SimVariablesModule, only: idm_context
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_setptr
  use STLStackIntModule
  use STLVecIntModule
  use CharacterStringModule
  implicit none
  private

  ! constants for memory allocation
  integer(I4B), parameter :: MAX_SECTIONS_PER_SLN = 40
  integer(I4B), public, parameter :: LEN_SECTION_TITLE = 128

  ! data structure to store measurements for a section
  type, private :: MeasuredSectionType
    character(len=LEN_SECTION_TITLE) :: title !< title to identify timed section in log
    real(DP) :: walltime !< walltime spent in section
    integer(I4B) :: count !< number of times section was entered
    integer(I4B) :: status !< =1 means section timer started, =0 otherwise
    integer(I4B) :: parent_id !< id of parent, or 0 when root
    type(STLVecInt) :: children !< ids of children
  end type MeasuredSectionType

  !> @brief A public type for profiling performance in the application.
  !! The ProfilerType is used to measure and record the performance of various
  !! parts of the application. It provides mechanisms to start, stop, and
  !< report on the performance metrics collected during execution.
  type, public :: ProfilerType
    real(DP) :: sim_start_time !< the simulation start time which lies before initialization of the profiler
    ! handles for the global simulation structure (with no simulation objects to store them)
    integer(I4B) :: tmr_run !< handle to timed section "Run"
    integer(I4B) :: tmr_init !< handle to timed section "Initialize"
    integer(I4B) :: tmr_update !< handle to timed section "Update"
    integer(I4B) :: tmr_finalize !< handle to timed section "Finalize"
    integer(I4B) :: tmr_prep_tstp !< handle to timed section "Prepare time step"
    integer(I4B) :: tmr_do_tstp !< handle to timed section "Do time step"
    integer(I4B) :: tmr_final_tstp !< handle to timed section "Finalize time step"
    integer(I4B) :: tmr_output !< handle to timed section "Write output"
    integer(I4B) :: tmr_nc_export !< handle to timed section "NetCDF export"
    ! private
    integer(I4B), private :: iout !< output unit number, typically simulation listing file
    integer(I4B), private :: pr_option !< 0 = NONE, 1 = SUMMARY, 2 = DETAIL
    integer(I4B), private :: nr_sections !< number of sections
    integer(I4B), private, dimension(3) :: top_three !< top three leaf sections based on walltime
    integer(I4B), private :: max_title_len !< maximum title length
    integer(I4B), private :: root_id !< currently only one root section is supported, this is the id
    type(MeasuredSectionType), dimension(:), pointer :: all_sections => null() !< all timed sections (up to MAX_NR_TIMED_SECTIONS)
    type(STLStackInt) :: callstack !< call stack of section ids
  contains
    procedure :: pre_init
    procedure :: initialize
    procedure :: start
    procedure :: stop
    procedure :: print
    procedure :: destroy
    ! private
    procedure, private :: set_print_option
    procedure, private :: add_section
    procedure, private :: print_section
    procedure, private :: print_total
    procedure, private :: aggregate_walltime
    procedure, private :: aggregate_counts
    procedure, private :: largest_title_length
    procedure, private :: sort_by_walltime
  end type ProfilerType

  type(ProfilerType), public :: g_prof !< the global timer object (to reduce trivial lines of code)

contains

  !< @brief To save the start time before initialization
  subroutine pre_init(this)
    class(ProfilerType) :: this

    call cpu_time(this%sim_start_time)

  end subroutine pre_init

  !< @brief Initialize the CPU timer object
  !<
  subroutine initialize(this)
    class(ProfilerType) :: this
    ! local
    character(len=LENMEMPATH) :: input_mempath
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: slntype
    character(len=:), pointer :: prprof
    integer(I4B) :: i
    integer(I4B) :: max_sections
    integer(I4B) :: nr_solutions

    this%tmr_run = -1
    this%tmr_init = -1
    this%tmr_update = -1
    this%tmr_finalize = -1
    this%tmr_prep_tstp = -1
    this%tmr_do_tstp = -1
    this%tmr_final_tstp = -1
    this%tmr_output = -1
    this%tmr_nc_export = -1

    call this%callstack%init()

    ! get nr of solutions from input context
    input_mempath = create_mem_path('SIM', 'NAM', idm_context)
    call mem_setptr(slntype, 'SLNTYPE', input_mempath)
    nr_solutions = size(slntype)
    call mem_setptr(prprof, 'PRPROF', input_mempath)
    call this%set_print_option(prprof)

    max_sections = MAX_SECTIONS_PER_SLN * nr_solutions
    allocate (this%all_sections(max_sections))
    do i = 1, max_sections
      this%all_sections(i)%title = "undefined"
      this%all_sections(i)%status = 0
      this%all_sections(i)%walltime = DZERO
      this%all_sections(i)%count = 0
      this%all_sections(i)%parent_id = 0
      call this%all_sections(i)%children%init()
    end do

    this%nr_sections = 0
    this%root_id = 0
    this%top_three = [0, 0, 0]

    ! start root section here with previously recorded start time
    if (this%pr_option > 0) then
      call this%start("Run", this%tmr_run)
      this%all_sections(this%tmr_run)%walltime = -this%sim_start_time
    end if

  end subroutine initialize

  !> @brief Add a new timed section to the tree,
  !! passing the parent id will add it as a child
  !< in the tree
  function add_section(this, title, parent_id) result(section_id)
    use SimModule, only: ustop
    class(ProfilerType) :: this
    character(len=*) :: title
    integer(I4B) :: parent_id
    integer(I4B) :: section_id

    ! increment to new section id
    this%nr_sections = this%nr_sections + 1
    section_id = this%nr_sections
    if (section_id > size(this%all_sections)) then
      write (*, *) "Internal error: Too many profiled sections, "&
        &"disable profiling to circumvent."
      call ustop()
    end if

    ! initialize new section
    this%all_sections(section_id)%title = title
    this%all_sections(section_id)%walltime = DZERO
    this%all_sections(section_id)%status = 0

    ! if parent, otherwise root section
    if (parent_id > 0) then
      ! add child to parent
      this%all_sections(section_id)%parent_id = parent_id
      call this%all_sections(parent_id)%children%push_back(section_id)
    else
      ! this is the root, assume there's only one!
      this%all_sections(section_id)%parent_id = 0
      this%root_id = section_id
    end if

  end function add_section

  !> @brief Start section timing, add when not exist yet (i.e. when id < 1)
  !<
  subroutine start(this, title, section_id)
    class(ProfilerType) :: this
    character(len=*) :: title
    integer(I4B) :: section_id
    ! local
    integer(I4B) :: parent_id
    real(DP) :: start_time
    type(MeasuredSectionType), pointer :: section

    if (this%pr_option == 0) return

    call cpu_time(start_time)

    if (section_id == -1) then
      ! add section if not exist
      parent_id = 0 ! root
      if (this%callstack%size() > 0) then
        parent_id = this%callstack%top()
      end if
      section_id = this%add_section(title, parent_id)
    end if
    call this%callstack%push(section_id)

    section => this%all_sections(section_id)
    section%count = section%count + 1
    section%status = 1
    section%walltime = section%walltime - start_time

  end subroutine start

  subroutine stop(this, section_id)
    class(ProfilerType) :: this
    integer(I4B) :: section_id
    ! local
    real(DP) :: end_time
    type(MeasuredSectionType), pointer :: section

    if (this%pr_option == 0) return

    call cpu_time(end_time)

    ! nett result (c.f. start(...)) is adding (dt = end_time - start_time)
    section => this%all_sections(section_id)
    section%status = 0
    section%walltime = section%walltime + end_time

    ! pop from call stack
    call this%callstack%pop()

  end subroutine stop

  subroutine print(this, output_unit)
    class(ProfilerType) :: this
    integer(I4B), intent(in) :: output_unit
    ! local
    integer(I4B) :: level, i, top_idx
    integer(I4B), dimension(:), allocatable :: sorted_idxs

    if (this%pr_option == 0) return

    this%iout = output_unit

    ! get top three leaf sections based on walltime
    top_idx = 1
    sorted_idxs = (/(i, i=1, this%nr_sections)/)
    call this%sort_by_walltime(sorted_idxs)
    do i = 1, this%nr_sections
      if (this%all_sections(sorted_idxs(i))%children%size == 0) then ! leaf node
        if (top_idx > 3) exit
        this%top_three(top_idx) = sorted_idxs(i)
        top_idx = top_idx + 1
      end if
    end do

    this%max_title_len = this%largest_title_length()

    if (this%pr_option > 1) then
      ! print timing call stack
      level = 0
      write (this%iout, '(/1x,a/)') &
        repeat('-', 18)//" Profiler: Call Stack "//repeat('-', 18)
      call this%print_section(this%root_id, level)
    end if

    ! print walltime per category from substring (if exist)
    ! note: the sections containing the substring should not be nested,
    !       otherwise the walltime will be counted multiple times
    write (this%iout, '(1x,a/)') &
      repeat('-', 20)//" Profiler: Totals "//repeat('-', 20)
    call this%print_total("Formulate")
    call this%print_total("Linear solve")
    call this%print_total("Calculate flows")
    call this%print_total("Calculate budgets")
    call this%print_total("Write output")
    call this%print_total("Parallel Solution")
    call this%print_total("MPI_WaitAll")
    write (this%iout, '(/1x,a/)') &
      repeat('-', 22)//" End Profiler "//repeat('-', 22)

  end subroutine print

  recursive subroutine print_section(this, section_id, level)
    use ArrayHandlersModule, only: ifind
    class(ProfilerType) :: this
    integer(I4B) :: section_id
    integer(I4B) :: level
    ! local
    integer(I4B) :: i, new_level, nr_padding, idx_top
    real(DP) :: percent
    type(MeasuredSectionType), pointer :: section
    character(len=:), allocatable :: title_padded
    character(len=LINELENGTH) :: top_marker

    section => this%all_sections(section_id)

    ! calculate percentage
    percent = 1.0_DP
    if (section%parent_id /= 0) then
      percent = section%walltime / this%all_sections(this%root_id)%walltime
    end if
    percent = percent * 100.0_DP

    ! determine if section should be marked as top three
    top_marker = ""
    idx_top = ifind(this%top_three, section_id)
    if (idx_top > 0) then
      nr_padding = max(0, 32 - level * 4)
      write (top_marker, '(a,i0)') repeat(" ", nr_padding)//"<== #", idx_top
    end if

    ! print section timing
    nr_padding = this%max_title_len - len_trim(section%title) + 2
    title_padded = trim(section%title)//":"//repeat(' ', nr_padding)
    write (this%iout, '(3a,f6.2,2a,f14.6,2a,i0,a,a)') " ", &
      repeat('....', level), "[", percent, "%] ", title_padded, &
      section%walltime, "s", " (", section%count, "x)", trim(top_marker)

    ! print children
    new_level = level + 1
    do i = 1, section%children%size
      call this%print_section(section%children%at(i), new_level)
    end do

    if (level == 0) write (this%iout, *)

  end subroutine print_section

  subroutine print_total(this, subtitle)
    class(ProfilerType) :: this
    character(len=*) :: subtitle
    ! local
    integer(I4B) :: count
    real(DP) :: walltime, percent
    integer(I4B) :: nr_padding
    character(len=:), allocatable :: title_padded

    ! get maximum length of title
    nr_padding = this%max_title_len - len_trim(subtitle)
    title_padded = trim(subtitle)//repeat(' ', nr_padding)

    count = this%aggregate_counts(subtitle)
    if (count > 0) then
      walltime = aggregate_walltime(this, subtitle)
      percent = (walltime / this%all_sections(this%root_id)%walltime) * 100.0_DP
      write (this%iout, '(2a,f6.2,3a,f14.6,2a,i0,a)') " ", "[", percent, &
        "%] ", title_padded, ": ", walltime, "s", " (", count, "x)"
    end if

  end subroutine print_total

  !> @brief Aggregate walltime over sections with a certain title
  !<
  function aggregate_walltime(this, title) result(walltime)
    class(ProfilerType) :: this
    character(len=*) :: title
    real(DP) :: walltime
    ! local
    integer(I4B) :: i

    walltime = DZERO
    do i = 1, this%nr_sections
      if (index(this%all_sections(i)%title, trim(title)) > 0) then
        walltime = walltime + this%all_sections(i)%walltime
      end if
    end do

  end function aggregate_walltime

  !> @brief Aggregate counts over sections with a certain title
  !<
  function aggregate_counts(this, title) result(counts)
    class(ProfilerType) :: this
    character(len=*) :: title
    integer(I4B) :: counts
    ! local
    integer(I4B) :: i

    counts = 0
    do i = 1, this%nr_sections
      if (index(this%all_sections(i)%title, trim(title)) > 0) then
        counts = counts + this%all_sections(i)%count
      end if
    end do

  end function aggregate_counts

  !> @brief Set the profile option from the user input
  !<
  subroutine set_print_option(this, profile_option)
    class(ProfilerType) :: this
    character(len=*), intent(in) :: profile_option

    select case (trim(profile_option))
    case ("NONE")
      this%pr_option = 0
    case ("SUMMARY")
      this%pr_option = 1
    case ("DETAIL")
      this%pr_option = 2
    case default
      this%pr_option = 0
    end select

  end subroutine set_print_option

  !> @brief Clean up the CPU timer object
  !<
  subroutine destroy(this)
    class(ProfilerType) :: this
    ! local
    integer(I4B) :: i

    call this%callstack%destroy()

    do i = 1, size(this%all_sections)
      call this%all_sections(i)%children%destroy()
    end do
    deallocate (this%all_sections)
    nullify (this%all_sections)

  end subroutine destroy

  !> @brief Calculate the largest title length
  !<
  function largest_title_length(this) result(max_length)
    class(ProfilerType) :: this
    integer(I4B) :: max_length
    integer(I4B) :: i

    max_length = 0
    do i = 1, this%nr_sections
      max_length = max(max_length, len_trim(this%all_sections(i)%title))
    end do

  end function largest_title_length

  !> @brief Sort section indexes based on walltime
  !<
  subroutine sort_by_walltime(this, idxs)
    class(ProfilerType) :: this
    integer(I4B), dimension(:), allocatable :: idxs !< array with unsorted section idxs
    integer(I4B) :: i, j, temp

    ! Simple bubble sort for demonstration purposes
    do i = 1, size(idxs) - 1
      do j = 1, size(idxs) - i
        if (this%all_sections(idxs(j))%walltime < &
            this%all_sections(idxs(j + 1))%walltime) then
          temp = idxs(j)
          idxs(j) = idxs(j + 1)
          idxs(j + 1) = temp
        end if
      end do
    end do

  end subroutine sort_by_walltime

end module ProfilerModule
