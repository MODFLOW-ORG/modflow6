!> @brief This module contains common time-varying property functionality
!!
!! This module contains methods implementing functionality common to both
!! time-varying hydraulic conductivity (TVK) and time-varying storage (TVS)
!! packages.
!!
!<
module TvBaseModule
  use BaseDisModule, only: DisBaseType
  use ConstantsModule, only: LINELENGTH, MAXCHARLEN, DZERO
  use GeomUtilModule, only: get_node
  use KindModule, only: I4B, DP, LGP
  use NumericalPackageModule, only: NumericalPackageType
  use SimModule, only: count_errors, store_error, store_error_filename
  use SimVariablesModule, only: errmsg
  use TdisModule, only: kper, nper, kstp
  use MemoryManagerModule, only: mem_setptr, get_isize

  implicit none

  private

  public :: TvBaseType
  public :: tvbase_da

  type, abstract, extends(NumericalPackageType) :: TvBaseType
    integer(I4B), dimension(:, :), pointer, contiguous :: cellid => null() !< input cellids of time varying value
    logical(LGP) :: ts_active = .false. !< is timeseries active in package
  contains
    procedure :: init
    procedure :: ar
    procedure :: rp
    procedure :: ad
    procedure :: da => tvbase_da
    procedure, private :: tvbase_allocate_scalars
    procedure, private :: tv_get_node
    procedure(ar_set_pointers), deferred :: ar_set_pointers
    procedure :: source_options => tvbase_source_options
    procedure :: source_package_options => tvbase_source_package_options
    procedure(apply_row_changes), deferred :: apply_row_changes
    procedure(set_changed_at), deferred :: set_changed_at
    procedure(reset_change_flags), deferred :: reset_change_flags
    procedure(validate_change), deferred :: validate_change
  end type TvBaseType

  abstract interface

    !> @brief Announce package and set pointers to variables
    !!
    !! Deferred procedure called by the TvBaseType code to announce the
    !! specific package version and set any required array and variable
    !! pointers from other packages.
    !<
    subroutine ar_set_pointers(this)
      ! -- modules
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
    end subroutine

    !> @brief Apply input column changes for period-data row n to node.
    !<
    subroutine apply_row_changes(this, n, node)
      ! -- modules
      use KindModule, only: I4B
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      integer(I4B), intent(in) :: n !< row index in period data arrays
      integer(I4B), intent(in) :: node !< reduced node number
    end subroutine

    !> @brief Mark property changes as having occurred at (kper, kstp)
    !!
    !! Deferred procedure called by the TvBaseType code when a property value
    !! change occurs at (kper, kstp).
    !<
    subroutine set_changed_at(this, kper, kstp)
      ! -- modules
      use KindModule, only: I4B
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      integer(I4B), intent(in) :: kper
      integer(I4B), intent(in) :: kstp
    end subroutine

    !> @brief Clear all per-node change flags
    !!
    !! Deferred procedure called by the TvBaseType code when a new time step
    !! commences, indicating that any previously set per-node property value
    !! change flags should be reset.
    !<
    subroutine reset_change_flags(this)
      ! -- modules
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
    end subroutine

    !> @brief Check that a given property value is valid
    !!
    !! Deferred procedure called by each derived type's apply_row_changes
    !! implementation after a property value change occurs. Performs any
    !! required validity checks on the value of the given variable at the
    !! given node. Perform any required updates to the property value if it
    !! is valid, or log an error if not.
    !<
    subroutine validate_change(this, n, varName)
      ! -- modules
      use KindModule, only: I4B
      import TvBaseType
      ! -- dummy
      class(TvBaseType) :: this
      integer(I4B), intent(in) :: n
      character(len=*), intent(in) :: varName
    end subroutine

  end interface

contains

  !> @brief Initialize the TvBaseType object
  !!
  !! Allocate and initialize data members of the object.
  !<
  subroutine init(this, name_model, pakname, ftype, mempath, inunit, iout)
    ! -- dummy
    class(TvBaseType) :: this
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: ftype
    character(len=*), intent(in) :: mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    call this%set_names(1, name_model, pakname, ftype, mempath)
    call this%tvbase_allocate_scalars()
    this%inunit = inunit
    this%iout = iout
  end subroutine init

  !> @brief Allocate scalar variables
  !!
  !! Allocate scalar data members of the object.
  !<
  subroutine tvbase_allocate_scalars(this)
    ! -- dummy
    class(TvBaseType) :: this
    !
    ! -- Call standard NumericalPackageType allocate scalars
    call this%NumericalPackageType%allocate_scalars()
  end subroutine tvbase_allocate_scalars

  !> @brief Source common options from the input memory path.
  !!
  !! Source common options and call derived package routine to source
  !! and log any package-specific options within the same block.
  !<
  subroutine tvbase_source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(TvBaseType) :: this
    ! -- locals
    integer(I4B) :: isize
    logical(LGP) :: found_print_input
    !
    write (this%iout, '(1x,a)') &
      'PROCESSING '//trim(adjustl(this%packName))//' OPTIONS'
    !
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, &
                       found_print_input)
    !
    if (found_print_input) then
      write (this%iout, '(4x,a)') 'TIME-VARYING INPUT WILL BE PRINTED.'
    end if
    !
    call get_isize('TS6_FILENAME', this%input_mempath, isize)
    if (isize > 0) this%ts_active = .true.
    !
    ! -- source package-specific options
    call this%source_package_options()
    !
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%packName))//' OPTIONS'
  end subroutine tvbase_source_options

  !> @brief Source package-specific options from the input memory path.
  !!
  !! Override in derived package to source and log package-specific options.
  !! Default implementation is a no-op.
  !<
  subroutine tvbase_source_package_options(this)
    ! -- dummy
    class(TvBaseType) :: this
    !
    ! -- no package-specific options in the base class
  end subroutine tvbase_source_package_options

  !> @brief Allocate and read static data for the package.
  !<
  subroutine ar(this, dis)
    ! -- dummy
    class(TvBaseType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    !
    this%dis => dis
    call this%ar_set_pointers()
    !
    call this%source_options()
    !
    ! -- set input mempath pointers
    call mem_setptr(this%cellid, 'CELLID', this%input_mempath)
    !
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine ar

  !> @brief Read and prepare stress period data for the package.
  !<
  subroutine rp(this)
    ! -- dummy
    class(TvBaseType) :: this
    ! -- local variables
    integer(I4B), pointer :: iper, nbound
    integer(I4B) :: n, node
    character(len=LINELENGTH) :: cellstr
    !
    ! -- check last loaded input period
    call mem_setptr(iper, 'IPER', this%input_mempath)
    if (iper /= kper) return
    !
    ! -- When timeseries are active, ad applies values at every time step
    if (this%ts_active) return
    !
    call mem_setptr(nbound, 'NBOUND', this%input_mempath)
    !
    ! -- Reset per-node property change flags
    call this%reset_change_flags()
    !
    do n = 1, nbound
      node = this%tv_get_node(n)
      if (node < 1 .or. node > this%dis%nodes) then
        call this%dis%noder_to_string(node, cellstr)
        write (errmsg, '(a,2(1x,a))') &
          'CELLID', trim(cellstr), 'is not in the active model domain.'
        call store_error(errmsg)
        cycle
      end if
      !
      call this%apply_row_changes(n, node)
    end do
    !
    ! -- Record that changes were made at the current stress period / time step
    if (nbound > 0) then
      call this%set_changed_at(kper, kstp)
    end if
    !
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine rp

  !> @brief Apply advanced values at each time step.
  !<
  subroutine ad(this)
    ! -- dummy
    class(TvBaseType) :: this
    ! -- local variables
    integer(I4B), pointer :: nbound
    integer(I4B) :: n, node
    !
    ! -- no-op when timeseries aren't active.
    if (.not. this%ts_active) return
    !
    call mem_setptr(nbound, 'NBOUND', this%input_mempath)
    if (nbound > 0) then
      ! -- Record that changes were made at the current time step
      call this%set_changed_at(kper, kstp)
      ! -- Reset node change flags
      call this%reset_change_flags()
      ! -- Apply row changes
      do n = 1, nbound
        node = this%tv_get_node(n)
        if (node < 1 .or. node > this%dis%nodes) cycle
        call this%apply_row_changes(n, node)
      end do
    end if
    !
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine ad

  !> @brief Deallocate package memory
  !!
  !! Deallocate package scalars and arrays.
  !<
  subroutine tvbase_da(this)
    ! -- dummy
    class(TvBaseType) :: this
    !
    nullify (this%cellid)
    call this%NumericalPackageType%da()
  end subroutine tvbase_da

  !> @brief Return the reduced node number for CELLID row n; caller must bounds-check.
  !<
  function tv_get_node(this, n) result(node)
    ! -- dummy
    class(TvBaseType) :: this
    integer(I4B), intent(in) :: n !< row index in the period-data arrays
    ! -- return
    integer(I4B) :: node
    ! -- local
    integer(I4B) :: nodeu
    !
    if (this%dis%ndim == 1) then
      nodeu = this%cellid(1, n)
    elseif (this%dis%ndim == 2) then
      nodeu = get_node(this%cellid(1, n), 1, &
                       this%cellid(2, n), &
                       this%dis%mshape(1), 1, &
                       this%dis%mshape(2))
    else
      nodeu = get_node(this%cellid(1, n), &
                       this%cellid(2, n), &
                       this%cellid(3, n), &
                       this%dis%mshape(1), &
                       this%dis%mshape(2), &
                       this%dis%mshape(3))
    end if
    node = this%dis%get_nodenumber(nodeu, 1)
  end function tv_get_node

end module TvBaseModule
