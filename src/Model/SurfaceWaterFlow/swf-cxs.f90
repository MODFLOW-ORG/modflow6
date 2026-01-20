! The SwfCxsType package is assigned for a model
! and can be used to calculate wetted area, wetted
! perimeter, hydraulic radius, composite roughness, etc.
! even if the user doesn't specify a CXS Package.
module SwfCxsModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, DZERO, DTWOTHIRDS
  use MemoryHelperModule, only: create_mem_path
  use MemoryManagerModule, only: mem_allocate
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename, count_errors
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType

  implicit none
  private
  public :: SwfCxsType, cxs_cr

  type, extends(NumericalPackageType) :: SwfCxsType

    ! provided as input
    integer(I4B), pointer :: nsections => null() !< number of cross section
    integer(I4B), pointer :: npoints => null() !< total number of cross-section points
    integer(I4B), dimension(:), pointer, contiguous :: idcxs => null() !< cross section id number, size nsections
    integer(I4B), dimension(:), pointer, contiguous :: nxspoints => null() !< number of cross section points for section, size nsections
    real(DP), dimension(:), pointer, contiguous :: xfraction => null() !< cross-section relative x distance, of size npoints
    real(DP), dimension(:), pointer, contiguous :: height => null() !< cross-section heights, of size npoints
    real(DP), dimension(:), pointer, contiguous :: manfraction => null() !< cross-section roughness data, of size npoints

    ! calculated from input
    integer(I4B), dimension(:), pointer, contiguous :: iacross => null() !< pointers to cross-section data for each section, of size nsections + 1

  contains

    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: source_options
    procedure :: log_options
    procedure :: source_dimensions
    procedure :: log_dimensions
    procedure :: source_packagedata
    procedure :: log_packagedata
    procedure :: check_packagedata
    procedure :: source_crosssectiondata
    procedure :: log_crosssectiondata
    procedure :: cxs_da
    procedure :: get_cross_section_info
    procedure :: get_area
    procedure :: get_wetted_perimeter => cxs_wetted_perimeter
    procedure :: get_roughness
    procedure :: get_conveyance => cxs_conveyance
    procedure :: get_hydraulic_radius
    procedure :: get_wetted_top_width
    procedure :: get_maximum_top_width
    procedure :: write_cxs_table

  end type SwfCxsType

contains

  !> @brief create package
  !<
  subroutine cxs_cr(pobj, name_model, input_mempath, inunit, iout, dis)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    type(SwfCxsType), pointer :: pobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(inout) :: dis !< the pointer to the discretization
    ! -- locals
    logical(LGP) :: found_fname
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'CXS --  CROSS SECTION PACKAGE, VERSION 1, 5/24/2023', &
       &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! -- Create the object
    allocate (pobj)

    ! -- create name and memory path
    call pobj%set_names(1, name_model, 'CXS', 'CXS')

    ! -- Allocate scalars
    call pobj%allocate_scalars()

    ! -- Set variables
    pobj%input_mempath = input_mempath
    pobj%inunit = inunit
    pobj%iout = iout
    pobj%dis => dis

    ! -- set name of input file
    call mem_set_value(pobj%input_fname, 'INPUT_FNAME', pobj%input_mempath, &
                       found_fname)

    ! -- check if package is enabled
    if (inunit > 0) then

      ! -- Print a message identifying the package.
      write (iout, fmtheader) input_mempath

      ! -- source options
      call pobj%source_options()

      ! -- source dimensions
      call pobj%source_dimensions()

      ! -- allocate arrays
      call pobj%allocate_arrays()

      ! -- source dimensions
      call pobj%source_packagedata()

      ! -- source dimensions
      call pobj%source_crosssectiondata()

    end if
  end subroutine cxs_cr

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the package. The base model
  !! allocate scalars method is also called.
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    ! -- dummy
    class(SwfCxsType) :: this
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate scalars
    call mem_allocate(this%nsections, 'NSECTIONS', this%memoryPath)
    call mem_allocate(this%npoints, 'NPOINTS', this%memoryPath)

    ! -- initialize
    this%nsections = 0
    this%npoints = 0
  end subroutine allocate_scalars

  !> @brief Copy options from IDM into package
  !<
  subroutine source_options(this)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    ! -- dummy
    class(SwfCxsType) :: this
    ! -- locals
    character(len=LENMEMPATH) :: idmMemoryPath
    type(SwfCxsParamFoundType) :: found
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'CXS', idm_context)
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%iprpak, 'PRINT_INPUT', idmMemoryPath, &
                       found%iprpak)
    !
    ! -- log values to list file
    if (this%iout > 0) then
      call this%log_options(found)
    end if
  end subroutine source_options

  !> @brief Write user options to list file
  !<
  subroutine log_options(this, found)
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    class(SwfCxsType) :: this
    type(SwfCxsParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting CXS Options'

    if (found%iprpak) then
      write (this%iout, '(4x,a)') 'Package information will be printed.'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting CXS Options'

  end subroutine log_options

  !> @brief Copy options from IDM into package
  !<
  subroutine source_dimensions(this)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    ! -- dummy
    class(SwfCxsType) :: this
    ! -- locals
    character(len=LENMEMPATH) :: idmMemoryPath
    type(SwfCxsParamFoundType) :: found
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'CXS', idm_context)
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%nsections, 'NSECTIONS', idmMemoryPath, &
                       found%nsections)
    call mem_set_value(this%npoints, 'NPOINTS', idmMemoryPath, &
                       found%npoints)
    !
    ! -- ensure nsections was found
    if (.not. found%nsections) then
      write (errmsg, '(a)') 'Error in DIMENSIONS block: NSECTIONS not found.'
      call store_error(errmsg)
    end if
    !
    ! -- ensure npoints was found
    if (.not. found%npoints) then
      write (errmsg, '(a)') 'Error in DIMENSIONS block: NPOINTS not found.'
      call store_error(errmsg)
    end if
    !
    ! -- log values to list file
    if (this%iout > 0) then
      call this%log_dimensions(found)
    end if
  end subroutine source_dimensions

  !> @brief Write user options to list file
  !<
  subroutine log_dimensions(this, found)
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    class(SwfCxsType) :: this
    type(SwfCxsParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting CXS Dimensions'

    if (found%nsections) then
      write (this%iout, '(4x,a)') 'NSECTIONS set from input file.'
    end if

    if (found%npoints) then
      write (this%iout, '(4x,a)') 'NPOINTS set from input file.'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting CXS Dimensions'

  end subroutine log_dimensions

  !> @brief allocate memory for arrays
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(SwfCxsType) :: this
    ! -- locals
    integer(I4B) :: n
    !
    ! -- arrays allocation
    call mem_allocate(this%idcxs, this%nsections, &
                      'IDCXS', this%memoryPath)
    call mem_allocate(this%nxspoints, this%nsections, &
                      'NXSPOINTS', this%memoryPath)
    call mem_allocate(this%xfraction, this%npoints, &
                      'XFRACTION', this%memoryPath)
    call mem_allocate(this%height, this%npoints, &
                      'HEIGHT', this%memoryPath)
    call mem_allocate(this%manfraction, this%npoints, &
                      'MANFRACTION', this%memoryPath)
    call mem_allocate(this%iacross, this%nsections + 1, &
                      'IACROSS', this%memoryPath)

    ! -- initialization
    do n = 1, this%nsections
      this%idcxs(n) = 0
      this%nxspoints(n) = 0
    end do
    do n = 1, this%npoints
      this%xfraction(n) = DZERO
      this%height(n) = DZERO
      this%manfraction(n) = DZERO
    end do
    do n = 1, this%nsections + 1
      this%iacross(n) = 0
    end do
  end subroutine allocate_arrays

  !> @brief Copy options from IDM into package
  !<
  subroutine source_packagedata(this)
    ! modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    ! dummy
    class(SwfCxsType) :: this
    ! locals
    character(len=LENMEMPATH) :: idmMemoryPath
    type(SwfCxsParamFoundType) :: found

    ! set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'CXS', idm_context)

    ! update defaults with idm sourced values
    call mem_set_value(this%idcxs, 'IDCXS', idmMemoryPath, &
                       found%idcxs)
    call mem_set_value(this%nxspoints, 'NXSPOINTS', idmMemoryPath, &
                       found%nxspoints)

    ! ensure idcxs was found
    if (.not. found%idcxs) then
      write (errmsg, '(a)') 'Error in PACKAGEDATA block: IDCXS not found.'
      call store_error(errmsg)
    end if

    ! ensure nxspoints was found
    if (.not. found%nxspoints) then
      write (errmsg, '(a)') 'Error in PACKAGEDATA block: NXSPOINTS not found.'
      call store_error(errmsg)
    end if

    ! log values to list file
    if (this%iout > 0) then
      call this%log_packagedata(found)
    end if

    ! Check to make sure package data is valid
    call this%check_packagedata()

    ! Calculate the iacross index array using nxspoints
    call calc_iacross(this%nxspoints, this%iacross)
  end subroutine source_packagedata

  !> @brief Calculate index pointer array iacross from nxspoints
  !<
  subroutine calc_iacross(nxspoints, iacross)
    integer(I4B), dimension(:), intent(in) :: nxspoints
    integer(I4B), dimension(:), intent(inout) :: iacross
    integer(I4B) :: n
    iacross(1) = 1
    do n = 1, size(nxspoints)
      iacross(n + 1) = iacross(n) + nxspoints(n)
    end do
  end subroutine calc_iacross

  !> @brief Check packagedata
  !<
  subroutine check_packagedata(this)
    ! dummy arguments
    class(SwfCxsType) :: this !< this instance
    ! local variables
    integer(I4B) :: i

    ! Check that all cross section IDs are in range
    do i = 1, size(this%idcxs)
      if (this%idcxs(i) <= 0 .or. this%idcxs(i) > this%nsections) then
        write (errmsg, '(a, i0, a)') &
        'IDCXS values must be greater than 0 and less than NSECTIONS.  &
        &Found ', this%idcxs(i), '.'
        call store_error(errmsg)
      end if
    end do

    ! Check that nxspoints are greater than one
    do i = 1, size(this%nxspoints)
      if (this%nxspoints(i) <= 1) then
        write (errmsg, '(a, i0, a, i0, a)') &
        'NXSPOINTS values must be greater than 1 for each cross section.  &
        &Found ', this%nxspoints(i), ' for cross section ', this%idcxs(i), '.'
        call store_error(errmsg)
      end if
    end do

    ! write summary of package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

  end subroutine check_packagedata

  !> @brief Write user packagedata to list file
  !<
  subroutine log_packagedata(this, found)
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    class(SwfCxsType) :: this
    type(SwfCxsParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting CXS Package Data'

    if (found%idcxs) then
      write (this%iout, '(4x,a)') 'IDCXS set from input file.'
    end if

    if (found%nxspoints) then
      write (this%iout, '(4x,a)') 'NXSPOINTS set from input file.'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting CXS Package Data'

  end subroutine log_packagedata

  !> @brief Copy options from IDM into package
  !<
  subroutine source_crosssectiondata(this)
    ! -- modules
    use KindModule, only: LGP
    use MemoryManagerExtModule, only: mem_set_value
    use SimVariablesModule, only: idm_context
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    ! -- dummy
    class(SwfCxsType) :: this
    ! -- locals
    character(len=LENMEMPATH) :: idmMemoryPath
    type(SwfCxsParamFoundType) :: found
    !
    ! -- set memory path
    idmMemoryPath = create_mem_path(this%name_model, 'CXS', idm_context)
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%xfraction, 'XFRACTION', idmMemoryPath, &
                       found%xfraction)
    call mem_set_value(this%height, 'HEIGHT', idmMemoryPath, &
                       found%height)
    call mem_set_value(this%manfraction, 'MANFRACTION', idmMemoryPath, &
                       found%manfraction)
    !
    ! -- ensure xfraction was found
    if (.not. found%xfraction) then
      write (errmsg, '(a)') &
        'Error in CROSSSECTIONDATA block: xfraction not found.'
      call store_error(errmsg)
    end if
    !
    ! -- ensure height was found
    if (.not. found%height) then
      write (errmsg, '(a)') &
        'Error in CROSSSECTIONDATA block: HEIGHT not found.'
      call store_error(errmsg)
    end if
    !
    ! -- ensure manfraction was found
    if (.not. found%manfraction) then
      write (errmsg, '(a)') &
        'Error in CROSSSECTIONDATA block: MANFRACTION not found.'
      call store_error(errmsg)
    end if
    !
    ! -- log values to list file
    if (this%iout > 0) then
      call this%log_crosssectiondata(found)
    end if
  end subroutine source_crosssectiondata

  !> @brief Write user packagedata to list file
  !<
  subroutine log_crosssectiondata(this, found)
    use SwfCxsInputModule, only: SwfCxsParamFoundType
    class(SwfCxsType) :: this
    type(SwfCxsParamFoundType), intent(in) :: found

    write (this%iout, '(1x,a)') 'Setting CXS Cross Section Data'

    if (found%xfraction) then
      write (this%iout, '(4x,a)') 'XFRACTION set from input file.'
    end if

    if (found%height) then
      write (this%iout, '(4x,a)') 'HEIGHT set from input file.'
    end if

    if (found%manfraction) then
      write (this%iout, '(4x,a)') 'MANFRACTION set from input file.'
    end if

    write (this%iout, '(1x,a,/)') 'End Setting CXS Cross Section Data'

  end subroutine log_crosssectiondata

  subroutine write_cxs_table(this, idcxs, width, slope, rough, unitconv)
    ! -- module
    use SortModule, only: qsort, unique_values
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs
    real(DP), intent(in) :: width
    real(DP), intent(in) :: slope
    real(DP), intent(in) :: rough
    real(DP), intent(in) :: unitconv
    ! -- local
    integer(I4B) :: ipt
    real(DP) :: d
    real(DP) :: a
    real(DP) :: rh
    real(DP) :: wp
    real(DP) :: r
    real(DP) :: c
    real(DP) :: q
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    real(DP), dimension(:), allocatable :: depths
    real(DP), dimension(:), allocatable :: depths_unique
    integer(I4B), dimension(:), allocatable :: indx

    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)

    if (npts > 0) then

      write (this%iout, *) 'Processing information for cross section ', idcxs
      write (this%iout, *) 'Depth Area WettedP HydRad Rough Conveyance Q'

      allocate (depths(npts))
      allocate (indx(size(depths)))

      depths(:) = this%height(:)
      call qsort(indx, depths)
      call unique_values(depths, depths_unique)

      do ipt = 1, size(depths_unique)
        d = depths_unique(ipt)
        a = this%get_area(idcxs, width, d)
        wp = this%get_wetted_perimeter(idcxs, width, d)
        rh = this%get_hydraulic_radius(idcxs, width, d, a)
        r = this%get_roughness(idcxs, width, d, rough, slope)
        c = this%get_conveyance(idcxs, width, d, rough)
        if (slope > DZERO) then
          q = unitconv * c * sqrt(slope)
        else
          q = DZERO
        end if
        write (this%iout, *) d, a, wp, rh, r, c, q
      end do

      deallocate (depths)
      deallocate (depths_unique)
      write (this%iout, *) 'Done processing information for cross section ', idcxs

    end if
  end subroutine write_cxs_table

  !> @brief deallocate memory
  !<
  subroutine cxs_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    use MemoryManagerExtModule, only: memorystore_remove
    use SimVariablesModule, only: idm_context
    ! -- dummy
    class(SwfCxsType) :: this
    !
    ! -- Deallocate input memory
    call memorystore_remove(this%name_model, 'CXS', idm_context)
    !
    ! -- Scalars
    call mem_deallocate(this%nsections)
    call mem_deallocate(this%npoints)
    !
    ! -- Deallocate arrays if the package was created
    !    from an input file
    if (this%inunit > 0) then
      call mem_deallocate(this%idcxs)
      call mem_deallocate(this%nxspoints)
      call mem_deallocate(this%xfraction)
      call mem_deallocate(this%height)
      call mem_deallocate(this%manfraction)
      call mem_deallocate(this%iacross)
    end if
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine cxs_da

  subroutine get_cross_section_info(this, idcxs, i0, i1, npts, icalcmeth)
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id number
    integer(I4B), intent(inout) :: i0 !< starting cross section point number
    integer(I4B), intent(inout) :: i1 !< ending cross section point number
    integer(I4B), intent(inout) :: npts !< number of points in cross section
    integer(I4B), intent(inout) :: icalcmeth !< calculation method for mannings roughness
    ! -- local
    !
    ! -- Return npts = 0 if this package does not have input file
    if (this%inunit == 0 .or. idcxs == 0) then
      npts = 0
      i0 = 1
      i1 = 1
      icalcmeth = 0
    else
      !
      ! -- If the cross section id is 0, then it is a hydraulically wide channel,
      !    and only width and rough are needed (not xfraction, height, and manfraction)
      if (idcxs > 0) then
        i0 = this%iacross(idcxs)
        i1 = this%iacross(idcxs + 1) - 1
      else
        i0 = 1
        i1 = 1
      end if
      ! set icalcmeth based on number of cross section points
      npts = i1 - i0 + 1
      icalcmeth = 0 ! linear composite mannings resistance
      if (npts > 4) then
        icalcmeth = 0 ! sum q by cross section segments
      end if
    end if
  end subroutine get_cross_section_info

  function get_area(this, idcxs, width, depth) result(area)
    ! -- modules
    use SwfCxsUtilsModule, only: get_cross_section_area
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    real(DP), intent(in) :: depth !< stage in reach
    ! -- local
    real(DP) :: area
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (npts == 0) then
      area = width * depth
    else
      area = get_cross_section_area(npts, &
                                    this%xfraction(i0:i1), &
                                    this%height(i0:i1), &
                                    width, depth)
    end if
  end function get_area

  function cxs_wetted_perimeter(this, idcxs, width, depth) result(wp)
    ! -- modules
    use SwfCxsUtilsModule, only: get_wetted_perimeter
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    real(DP), intent(in) :: depth !< stage in reach
    ! -- local
    real(DP) :: wp
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (npts == 0) then
      wp = width
    else
      wp = get_wetted_perimeter(npts, &
                                this%xfraction(i0:i1), &
                                this%height(i0:i1), &
                                width, depth)
    end if
  end function cxs_wetted_perimeter

  function get_roughness(this, idcxs, width, depth, rough, &
                         slope) result(roughc)
    ! -- modules
    use SwfCxsUtilsModule, only: calc_composite_roughness
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    real(DP), intent(in) :: depth !< stage in reach
    real(DP), intent(in) :: rough !< mannings value provided for the reach
    real(DP), intent(in) :: slope !< slope value provided for the reach
    ! -- local
    real(DP) :: roughc !< calculated composite roughness
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (npts == 0) then
      roughc = rough
    else
      roughc = calc_composite_roughness(npts, &
                                        depth, &
                                        width, &
                                        rough, &
                                        slope, &
                                        this%xfraction(i0:i1), &
                                        this%height(i0:i1), &
                                        this%manfraction(i0:i1), &
                                        icalcmeth)
    end if
  end function get_roughness

  !> @brief Calculate and return conveyance
  !!
  !! Conveyance = area * hydraulic_radius ** (2/3) / mannings_roughness
  !! If idcxs = 0 (no cross section specified) then reach is
  !< hydraulically wide and hydraulic radius is equal to depth.
  function cxs_conveyance(this, idcxs, width, depth, &
                          rough) result(conveyance)
    ! -- modules
    use SwfCxsUtilsModule, only: get_conveyance
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    real(DP), intent(in) :: depth !< stage in reach
    real(DP), intent(in) :: rough !< mannings value provided for the reach
    ! -- return
    real(DP) :: conveyance !< calculated composite roughness
    ! -- local
    real(DP) :: a
    real(DP) :: rh
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (npts == 0) then
      a = depth * width
      rh = depth
      conveyance = a * rh**DTWOTHIRDS / rough
    else
      conveyance = get_conveyance(npts, &
                                  this%xfraction(i0:i1), &
                                  this%height(i0:i1), &
                                  this%manfraction(i0:i1), &
                                  width, rough, depth)
    end if
  end function cxs_conveyance

  function get_hydraulic_radius(this, idcxs, width, depth, area) result(r)
    ! -- modules
    use SwfCxsUtilsModule, only: get_hydraulic_radius_xf
    ! -- dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    real(DP), intent(in) :: depth !< stage in reach
    real(DP), intent(in), optional :: area !< area of the reach
    ! -- local
    real(DP) :: r !< calculated hydraulic radius
    real(DP) :: a
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (present(area)) then
      a = area
    else
      a = this%get_area(idcxs, width, depth)
    end if
    if (npts == 0) then
      r = a / width
    else
      r = get_hydraulic_radius_xf(npts, &
                                  this%xfraction(i0:i1), &
                                  this%height(i0:i1), &
                                  width, depth)
    end if
  end function get_hydraulic_radius

  function get_wetted_top_width(this, idcxs, width, depth) result(r)
    ! modules
    use SwfCxsUtilsModule, only: get_wetted_topwidth
    ! dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    real(DP), intent(in) :: depth !< stage in reach
    ! local
    real(DP) :: r !< calculated hydraulic radius
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (npts == 0) then
      r = width
    else
      r = get_wetted_topwidth(npts, this%xfraction(i0:i1), &
                              this%height(i0:i1), width, depth)
    end if
  end function get_wetted_top_width

  function get_maximum_top_width(this, idcxs, width) result(r)
    ! modules
    use SwfCxsUtilsModule, only: get_saturated_topwidth
    ! dummy
    class(SwfCxsType) :: this
    integer(I4B), intent(in) :: idcxs !< cross section id
    real(DP), intent(in) :: width !< width in reach
    ! local
    real(DP) :: r !< calculated hydraulic radius
    integer(I4B) :: i0
    integer(I4B) :: i1
    integer(I4B) :: npts
    integer(I4B) :: icalcmeth
    call this%get_cross_section_info(idcxs, i0, i1, npts, icalcmeth)
    if (npts == 0) then
      r = width
    else
      r = get_saturated_topwidth(npts, this%xfraction(i0:i1), width)
    end if
  end function get_maximum_top_width

end module SwfCxsModule
