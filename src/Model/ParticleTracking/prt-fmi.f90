module PrtFmiModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, LENAUXNAME, LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use FlowModelInterfaceModule, only: FlowModelInterfaceType
  use BaseDisModule, only: DisBaseType
  use BudgetObjectModule, only: BudgetObjectType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate

  implicit none
  private
  public :: PrtFmiType
  public :: fmi_cr
  public :: IFLOWFACE_TOP, IFLOWFACE_BOTTOM

  character(len=LENPACKAGENAME) :: text = '    PRTFMI'

  !> @brief IFLOWFACE numbers for top and bottom faces
  enum, bind(c)
    enumerator :: IFLOWFACE_TOP = -1
    enumerator :: IFLOWFACE_BOTTOM = -2
  end enum

  type, extends(FlowModelInterfaceType) :: PrtFmiType
    private
    integer(I4B), pointer, public :: max_faces => null() !< max number of 3d cell faces
    real(DP), dimension(:), pointer, contiguous, public :: &
      SourceFlows => null() !< cell source flows array
    real(DP), dimension(:), pointer, contiguous, public :: &
      SinkFlows => null() !< cell sink flows array
    real(DP), dimension(:), pointer, contiguous, public :: &
      StorageFlows => null() !< cell storage flows array
    real(DP), dimension(:, :), pointer, contiguous, public :: &
      BoundaryFlows => null() !< cell boundary flows array
    integer(I4B), dimension(:), pointer, contiguous, public :: &
      BoundaryFaces => null() !< bitmask of assigned boundary faces
    ! -- Node/connection maps to GWF's reduced numbering. Unassociated
    !    unless PRT's active domain is a strict subset of GWF's
    integer(I4B), dimension(:), pointer, contiguous, public :: &
      noder_gwf2prt => null() !< maps GWF noder to PRT noder (0 if PRT-inactive)
    integer(I4B), dimension(:), pointer, contiguous, public :: &
      noder_prt2gwf => null() !< maps PRT noder to GWF noder
    integer(I4B), dimension(:), pointer, contiguous, public :: &
      ipos_prt2gwf => null() !< maps PRT ia/ja position to GWF ia/ja position

  contains

    procedure :: fmi_ad
    procedure :: fmi_df => prtfmi_df
    procedure :: fmi_da => prtfmi_da
    procedure :: allocate_scalars => prtfmi_allocate_scalars
    procedure :: allocate_arrays => prtfmi_allocate_arrays
    procedure, private :: accumulate_flows
    procedure :: get_gwfflowja
    procedure :: mark_boundary_face
    procedure :: is_boundary_face
    procedure :: is_net_out_boundary_face
    procedure, private :: iflowface_to_icellface

  end type PrtFmiType

contains

  !> @brief Create a new PrtFmi object
  subroutine fmi_cr(fmiobj, name_model, input_mempath, inunit, iout)
    ! dummy
    type(PrtFmiType), pointer :: fmiobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(inout) :: inunit
    integer(I4B), intent(in) :: iout

    ! Create the object
    allocate (fmiobj)

    ! create name and memory path
    call fmiobj%set_names(1, name_model, 'FMI', 'FMI', input_mempath)
    fmiobj%text = text

    ! Allocate scalars
    call fmiobj%allocate_scalars()

    ! Set variables
    fmiobj%inunit = inunit
    fmiobj%iout = iout

    ! Assign dependent variable label
    fmiobj%depvartype = 'TRACKS          '

  end subroutine fmi_cr

  !> @brief Time step advance
  subroutine fmi_ad(this)
    ! modules
    use ConstantsModule, only: DHDRY
    ! dummy
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: n, ng
    character(len=15) :: nodestr
    character(len=*), parameter :: fmtdry = &
     &"(/1X,'WARNING: DRY CELL ENCOUNTERED AT ',a,';  RESET AS INACTIVE')"
    character(len=*), parameter :: fmtrewet = &
     &"(/1X,'DRY CELL REACTIVATED AT ', a)"

    ! Set flag to indicated that flows are being updated.  For the case where
    ! flows may be reused (only when flows are read from a file) then set
    ! the flag to zero to indicated that flows were not updated
    this%iflowsupdated = 1

    ! If reading flows from a budget file, read the next set of records
    if (this%iubud /= 0) call this%advance_bfr()

    ! If reading heads from a head file, read the next set of records
    if (this%iuhds /= 0) call this%advance_hfr()

    ! If mover flows are being read from file, read the next set of records
    if (this%iumvr /= 0) &
      call this%mvrbudobj%bfr_advance(this%dis, this%iout)

    ! Accumulate flows
    call this%accumulate_flows()

    ! if flow cell is dry, then set this%ibound = 0
    do n = 1, this%dis%nodes
      ! noder_prt2gwf translates to GWF numbering if GWF and PRT
      ! have different active domains
      if (associated(this%noder_prt2gwf)) then
        ng = this%noder_prt2gwf(n)
      else
        ng = n
      end if

      ! Calculate the ibound-like array that has 0 if saturation
      ! is zero and 1 otherwise
      if (this%gwfsat(ng) > DZERO) then
        this%ibdgwfsat0(n) = 1
      else
        this%ibdgwfsat0(n) = 0
      end if

      ! Check if active model cell is inactive for flow
      if (this%ibound(n) > 0) then
        if (this%gwfhead(ng) == DHDRY) then
          ! cell should be made inactive
          this%ibound(n) = 0
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtdry) trim(nodestr)
        end if
      end if

      ! Convert dry model cell to active if flow has rewet
      if (this%ibound(n) == 0) then
        if (this%gwfhead(ng) /= DHDRY) then
          ! cell is now wet
          this%ibound(n) = 1
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtrewet) trim(nodestr)
        end if
      end if
    end do

  end subroutine fmi_ad

  !> @brief Define the flow model interface
  subroutine prtfmi_df(this, dis, idryinactive)
    class(PrtFmiType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: idryinactive

    call this%FlowModelInterfaceType%fmi_df(dis, idryinactive)

    this%max_faces = this%dis%get_max_npolyverts() + 2
    if (this%max_faces > 32) then
      write (errmsg, '(a,i0,a,i0,a)') &
        'DISV grid contains a cell with ', this%max_faces - 2, &
        ' lateral faces. Cells may have at most 30 lateral faces.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
      return
    end if

  end subroutine prtfmi_df

  !> @brief Allocate scalars
  subroutine prtfmi_allocate_scalars(this)
    class(PrtFmiType) :: this

    call this%FlowModelInterfaceType%allocate_scalars()

    call mem_allocate(this%max_faces, 'MAX_FACES', this%memoryPath)
    this%max_faces = 0

  end subroutine prtfmi_allocate_scalars

  !> @brief Allocate arrays
  subroutine prtfmi_allocate_arrays(this, nodes)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: nodes

    ! allocate parent arrays
    call this%FlowModelInterfaceType%allocate_arrays(nodes)

    call mem_allocate(this%StorageFlows, nodes, &
                      'STORAGEFLOWS', this%memoryPath)
    call mem_allocate(this%SourceFlows, nodes, &
                      'SOURCEFLOWS', this%memoryPath)
    call mem_allocate(this%SinkFlows, nodes, &
                      'SINKFLOWS', this%memoryPath)
    call mem_allocate(this%BoundaryFlows, nodes, this%max_faces, &
                      'BOUNDARYFLOWS', this%memoryPath)
    call mem_allocate(this%BoundaryFaces, nodes, &
                      'BOUNDARYFACES', this%memoryPath)

  end subroutine prtfmi_allocate_arrays

  !> @brief Deallocate memory
  subroutine prtfmi_da(this)
    class(PrtFmiType) :: this

    call mem_deallocate(this%max_faces)
    call mem_deallocate(this%StorageFlows)
    call mem_deallocate(this%SourceFlows)
    call mem_deallocate(this%SinkFlows)
    call mem_deallocate(this%BoundaryFlows)
    call mem_deallocate(this%BoundaryFaces)
    if (associated(this%noder_gwf2prt)) call mem_deallocate(this%noder_gwf2prt)
    if (associated(this%noder_prt2gwf)) call mem_deallocate(this%noder_prt2gwf)
    if (associated(this%ipos_prt2gwf)) call mem_deallocate(this%ipos_prt2gwf)

    call this%FlowModelInterfaceType%fmi_da()

  end subroutine prtfmi_da

  !> @brief Accumulate flows
  subroutine accumulate_flows(this)
    ! dummy
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: j, i, ip, ib
    integer(I4B) :: iflowface, iauxiflowface, icellface
    real(DP) :: qbnd
    character(len=LENAUXNAME) :: auxname
    integer(I4B) :: naux

    this%StorageFlows = DZERO
    if (associated(this%noder_prt2gwf)) then
      ! PRT's active domain is a strict subset of GWF's: gwfstrgss/gwfstrgsy
      ! are sized to GWF's node numbering, so gather through the node map
      ! (a vector subscript, still a whole-array operation)
      if (this%igwfstrgss /= 0) &
        this%StorageFlows = this%StorageFlows + &
                            this%gwfstrgss(this%noder_prt2gwf)
      if (this%igwfstrgsy /= 0) &
        this%StorageFlows = this%StorageFlows + &
                            this%gwfstrgsy(this%noder_prt2gwf)
    else
      if (this%igwfstrgss /= 0) &
        this%StorageFlows = this%StorageFlows + this%gwfstrgss
      if (this%igwfstrgsy /= 0) &
        this%StorageFlows = this%StorageFlows + this%gwfstrgsy
    end if

    this%SourceFlows = DZERO
    this%SinkFlows = DZERO
    this%BoundaryFlows = DZERO
    this%BoundaryFaces = 0
    do ip = 1, this%nflowpack
      iauxiflowface = 0
      naux = this%gwfpackages(ip)%naux
      if (naux > 0) then
        do j = 1, naux
          auxname = this%gwfpackages(ip)%auxname(j)
          if (trim(adjustl(auxname)) == "IFLOWFACE") then
            iauxiflowface = j
            exit
          end if
        end do
      end if
      do ib = 1, this%gwfpackages(ip)%nbound
        i = this%gwfpackages(ip)%nodelist(ib)
        ! translate to PRT node numbering if needed
        if (associated(this%noder_gwf2prt)) i = this%noder_gwf2prt(i)
        if (i <= 0) cycle
        if (this%ibound(i) <= 0) cycle
        qbnd = this%gwfpackages(ip)%get_flow(ib)
        ! todo, after initial release: default iflowface values for different packages
        iflowface = 0
        icellface = 0
        if (iauxiflowface > 0) then
          iflowface = NINT(this%gwfpackages(ip)%auxvar(iauxiflowface, ib))
          icellface = this%iflowface_to_icellface(iflowface)
        end if
        if (icellface > 0) then
          call this%mark_boundary_face(i, icellface)
          this%BoundaryFlows(i, icellface) = &
            this%BoundaryFlows(i, icellface) + qbnd
        else if (qbnd .gt. DZERO) then
          this%SourceFlows(i) = this%SourceFlows(i) + qbnd
        else if (qbnd .lt. DZERO) then
          this%SinkFlows(i) = this%SinkFlows(i) + qbnd
        end if
      end do
    end do

  end subroutine accumulate_flows

  !> @brief Get the flow across connection ipos, given in this model's own
  !! reduced ia/ja numbering, translating through GWF's connection numbering
  !! first if PRT's active domain is a strict subset of GWF's.
  function get_gwfflowja(this, ipos) result(q)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: ipos !< connection position in PRT's dis%con
    real(DP) :: q
    ! local
    integer(I4B) :: ipos_gwf

    ipos_gwf = ipos
    if (associated(this%ipos_prt2gwf)) ipos_gwf = this%ipos_prt2gwf(ipos)
    q = this%gwfflowja(ipos_gwf)
  end function get_gwfflowja

  !> @brief Mark a face as a boundary face.
  subroutine mark_boundary_face(this, ic, icellface)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: icellface !< cell face number
    ! local
    integer(I4B) :: bit_pos

    if (ic <= 0 .or. ic > this%dis%nodes) then
      print *, 'Invalid cell number: ', ic
      print *, 'Expected a value in range [1, ', this%dis%nodes, ']'
      call pstop(1)
    end if
    if (icellface <= 0) then
      print *, 'Invalid face number: ', icellface
      print *, 'Expected a value in range [1, ', this%max_faces, ']'
      call pstop(1)
    end if
    bit_pos = icellface - 1 ! bit position 0-based
    if (bit_pos < 0 .or. bit_pos > 31) then
      print *, 'Invalid bitmask position: ', bit_pos
      print *, 'Expected a value in range [0, 31]'
      call pstop(1)
    end if
    this%BoundaryFaces(ic) = ibset(this%BoundaryFaces(ic), bit_pos)
  end subroutine mark_boundary_face

  !> @brief Check if a face is assigned to a boundary package.
  function is_boundary_face(this, ic, icellface) result(is_boundary)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: icellface !< cell face number
    logical(LGP) :: is_boundary
    ! local
    integer(I4B) :: bit_pos

    is_boundary = .false.
    if (ic <= 0 .or. ic > this%dis%nodes) then
      print *, 'Invalid cell number: ', ic
      print *, 'Expected a value in range [1, ', this%dis%nodes, ']'
      call pstop(1)
    end if
    if (icellface <= 0) then
      print *, 'Invalid face number: ', icellface
      print *, 'Expected a value in range [1, ', this%max_faces, ']'
      call pstop(1)
    end if
    bit_pos = icellface - 1 ! bit position 0-based
    if (bit_pos < 0 .or. bit_pos > 31) then
      print *, 'Invalid bitmask position: ', bit_pos
      print *, 'Expected a value in range [0, 31]'
      call pstop(1)
    end if
    is_boundary = btest(this%BoundaryFaces(ic), bit_pos)
  end function is_boundary_face

  !> @brief Check if a face is an assigned boundary with net outflow.
  function is_net_out_boundary_face(this, ic, icellface) &
    result(is_net_out_boundary)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: icellface !< cell face number
    logical(LGP) :: is_net_out_boundary

    is_net_out_boundary = .false.
    if (.not. this%is_boundary_face(ic, icellface)) return
    if (this%BoundaryFlows(ic, icellface) < DZERO) &
      is_net_out_boundary = .true.
  end function is_net_out_boundary_face

  !> @brief Convert an iflowface number to a cell face number.
  !! Maps bottom (-2) -> max_faces - 1, top (-1) -> max_faces.
  function iflowface_to_icellface(this, iflowface) result(icellface)
    class(PrtFmiType), intent(inout) :: this
    integer(I4B), intent(in) :: iflowface
    integer(I4B) :: icellface

    icellface = iflowface
    if (icellface < 0) icellface = icellface + this%max_faces - IFLOWFACE_TOP
  end function iflowface_to_icellface

end module PrtFmiModule
