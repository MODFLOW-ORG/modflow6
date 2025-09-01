module PrtFmiModule

  use KindModule, only: DP, I4B, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, LENAUXNAME, LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use FlowModelInterfaceModule, only: FlowModelInterfaceType
  use BaseDisModule, only: DisBaseType
  use BudgetObjectModule, only: BudgetObjectType

  implicit none
  private
  public :: PrtFmiType
  public :: fmi_cr

  character(len=LENPACKAGENAME) :: text = '    PRTFMI'

  type, extends(FlowModelInterfaceType) :: PrtFmiType

    integer(I4B) :: max_faces !< maximum number of faces for grid cell polygons
    double precision, allocatable, public :: SourceFlows(:) ! cell source flows array
    double precision, allocatable, public :: SinkFlows(:) ! cell sink flows array
    double precision, allocatable, public :: StorageFlows(:) ! cell storage flows array
    double precision, allocatable, public :: BoundaryFlows(:) ! cell boundary flows array
    integer(I4B), allocatable, public :: BoundaryFaces(:) ! bitmask of assigned boundary faces

  contains

    procedure :: fmi_ad
    procedure :: fmi_df => prtfmi_df
    procedure, private :: accumulate_flows
    procedure :: mark_boundary_face
    procedure :: is_boundary_face

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
    !
    ! Create the object
    allocate (fmiobj)
    !
    ! create name and memory path
    call fmiobj%set_names(1, name_model, 'FMI', 'FMI', input_mempath)
    fmiobj%text = text
    !
    ! Allocate scalars
    call fmiobj%allocate_scalars()
    !
    ! Set variables
    fmiobj%inunit = inunit
    fmiobj%iout = iout
    !
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
    integer(I4B) :: n
    character(len=15) :: nodestr
    character(len=*), parameter :: fmtdry = &
     &"(/1X,'WARNING: DRY CELL ENCOUNTERED AT ',a,';  RESET AS INACTIVE')"
    character(len=*), parameter :: fmtrewet = &
     &"(/1X,'DRY CELL REACTIVATED AT ', a)"
    !
    ! Set flag to indicated that flows are being updated.  For the case where
    ! flows may be reused (only when flows are read from a file) then set
    ! the flag to zero to indicated that flows were not updated
    this%iflowsupdated = 1
    !
    ! If reading flows from a budget file, read the next set of records
    if (this%iubud /= 0) then
      call this%advance_bfr()
    end if
    !
    ! If reading heads from a head file, read the next set of records
    if (this%iuhds /= 0) then
      call this%advance_hfr()
    end if
    !
    ! If mover flows are being read from file, read the next set of records
    if (this%iumvr /= 0) then
      call this%mvrbudobj%bfr_advance(this%dis, this%iout)
    end if
    !
    ! Accumulate flows
    call this%accumulate_flows()
    !
    ! if flow cell is dry, then set this%ibound = 0
    do n = 1, this%dis%nodes
      !
      ! Calculate the ibound-like array that has 0 if saturation
      ! is zero and 1 otherwise
      if (this%gwfsat(n) > DZERO) then
        this%ibdgwfsat0(n) = 1
      else
        this%ibdgwfsat0(n) = 0
      end if
      !
      ! Check if active model cell is inactive for flow
      if (this%ibound(n) > 0) then
        if (this%gwfhead(n) == DHDRY) then
          ! cell should be made inactive
          this%ibound(n) = 0
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtdry) trim(nodestr)
        end if
      end if
      !
      ! Convert dry model cell to active if flow has rewet
      if (this%ibound(n) == 0) then
        if (this%gwfhead(n) /= DHDRY) then
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
    ! modules
    use SimModule, only: store_error
    ! dummy
    class(PrtFmiType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), intent(in) :: idryinactive
    !
    ! Call parent class define
    call this%FlowModelInterfaceType%fmi_df(dis, idryinactive)
    !
    ! Allocate arrays
    this%max_faces = this%dis%get_max_npolyverts() + 2
    allocate (this%StorageFlows(this%dis%nodes))
    allocate (this%SourceFlows(this%dis%nodes))
    allocate (this%SinkFlows(this%dis%nodes))
    allocate (this%BoundaryFlows(this%dis%nodes * this%max_faces))
    allocate (this%BoundaryFaces(this%dis%nodes))

  end subroutine prtfmi_df

  !> @brief Accumulate flows
  subroutine accumulate_flows(this)
    implicit none
    ! dummy
    class(PrtFmiType) :: this
    ! local
    integer(I4B) :: j, i, ip, ib
    integer(I4B) :: ioffset, iflowface, iauxiflowface, iface
    real(DP) :: qbnd
    character(len=LENAUXNAME) :: auxname
    integer(I4B) :: naux

    this%StorageFlows = DZERO
    if (this%igwfstrgss /= 0) &
      this%StorageFlows = this%StorageFlows + &
                          this%gwfstrgss
    if (this%igwfstrgsy /= 0) &
      this%StorageFlows = this%StorageFlows + &
                          this%gwfstrgsy

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
        if (i <= 0) cycle
        if (this%ibound(i) <= 0) cycle
        qbnd = this%gwfpackages(ip)%get_flow(ib)
        ! todo, after initial release: default iflowface values for different packages
        iflowface = 0 ! iflowface number
        iface = 0 ! internal face number
        if (iauxiflowface > 0) then
          iflowface = NINT(this%gwfpackages(ip)%auxvar(iauxiflowface, ib))
          iface = iflowface
          ! maps bot -2 -> max_faces - 1, top -1 -> max_faces
          if (iface < 0) iface = iface + this%max_faces + 1
        end if
        if (iface > 0) then
          call this%mark_boundary_face(i, iface)
          ioffset = (i - 1) * this%max_faces
          this%BoundaryFlows(ioffset + iface) = &
            this%BoundaryFlows(ioffset + iface) + qbnd
        else if (qbnd .gt. DZERO) then
          this%SourceFlows(i) = this%SourceFlows(i) + qbnd
        else if (qbnd .lt. DZERO) then
          this%SinkFlows(i) = this%SinkFlows(i) + qbnd
        end if
      end do
    end do

  end subroutine accumulate_flows

  !> @brief Mark a face as a boundary face.
  subroutine mark_boundary_face(this, ic, iface)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: iface !< face number
    ! local
    integer(I4B) :: bit_pos

    if (ic <= 0 .or. ic > size(this%BoundaryFaces)) return
    if (iface == 0) return
    bit_pos = iface - 1 ! bit position 0-based
    if (bit_pos < 0 .or. bit_pos > 31) then
      print *, 'Invalid bitmask position: ', iface
      print *, 'Expected a value in range [0, 31]'
      call pstop(1)
    end if
    this%BoundaryFaces(ic) = ibset(this%BoundaryFaces(ic), bit_pos)
  end subroutine mark_boundary_face

  !> @brief Check if a face is assigned to a boundary package.
  function is_boundary_face(this, ic, iface) result(is_boundary)
    class(PrtFmiType) :: this
    integer(I4B), intent(in) :: ic !< node number (reduced)
    integer(I4B), intent(in) :: iface !< face number
    logical(LGP) :: is_boundary
    ! local
    integer(I4B) :: bit_pos

    is_boundary = .false.
    if (ic <= 0 .or. ic > size(this%BoundaryFaces)) return
    if (iface == 0) return
    bit_pos = iface - 1 ! bit position 0-based
    if (bit_pos < 0 .or. bit_pos > 31) then
      print *, 'Invalid bitmask position: ', iface
      print *, 'Expected a value in range [0, 31]'
      call pstop(1)
    end if
    is_boundary = btest(this%BoundaryFaces(ic), bit_pos)
  end function is_boundary_face

end module PrtFmiModule
