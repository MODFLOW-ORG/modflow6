module TspAdvModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DONE, DZERO, DHALF, DTWO, DNODATA, DPREC, &
                             LINELENGTH
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use TspFmiModule, only: TspFmiType
  use TspAdvOptionsModule, only: TspAdvOptionsType
  use MatrixBaseModule

  implicit none
  private
  public :: TspAdvType
  public :: adv_cr

  type, extends(NumericalPackageType) :: TspAdvType

    integer(I4B), pointer :: iadvwt => null() !< advection scheme (0 up, 1 central, 2 tvd)
    real(DP), pointer :: ats_percel => null() !< user-specified fractional number of cells advection can move a particle during one time step
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy

  contains

    procedure :: adv_df
    procedure :: adv_ar
    procedure :: adv_dt
    procedure :: adv_fc
    procedure :: adv_cq
    procedure :: adv_da

    procedure :: allocate_scalars
    procedure, private :: read_options
    procedure, private :: advqtvd
    procedure, private :: advtvd_bd
    procedure, private :: advqtvd_experimental
    procedure, private :: compute_cell_gradient
    procedure :: adv_weight
    procedure :: advtvd

  end type TspAdvType

contains

  !> @ brief Create a new ADV object
  !!
  !!  Create a new ADV package
  !<
  subroutine adv_cr(advobj, name_model, inunit, iout, fmi, eqnsclfac)
    ! -- dummy
    type(TspAdvType), pointer :: advobj
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    type(TspFmiType), intent(in), target :: fmi
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    !
    ! -- Create the object
    allocate (advobj)
    !
    ! -- create name and memory path
    call advobj%set_names(1, name_model, 'ADV', 'ADV')
    !
    ! -- Allocate scalars
    call advobj%allocate_scalars()
    !
    ! -- Set variables
    advobj%inunit = inunit
    advobj%iout = iout
    advobj%fmi => fmi
    advobj%eqnsclfac => eqnsclfac
  end subroutine adv_cr

  !> @brief Define ADV object
  !!
  !! Define the ADV package
  !<
  subroutine adv_df(this, adv_options)
    ! -- dummy
    class(TspAdvType) :: this
    type(TspAdvOptionsType), optional, intent(in) :: adv_options !< the optional options, for when not constructing from file
    ! -- local
    character(len=*), parameter :: fmtadv = &
      "(1x,/1x,'ADV-- ADVECTION PACKAGE, VERSION 1, 8/25/2017', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Read or set advection options
    if (.not. present(adv_options)) then
      !
      ! -- Initialize block parser (adv has no define, so it's
      ! not done until here)
      call this%parser%Initialize(this%inunit, this%iout)
      !
      ! --print a message identifying the advection package.
      write (this%iout, fmtadv) this%inunit
      !
      ! --read options from file
      call this%read_options()
    else
      !
      ! --set options from input arg
      this%iadvwt = adv_options%iAdvScheme
    end if
  end subroutine adv_df

  !> @brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the ADV package.
  !<
  subroutine adv_ar(this, dis, ibound)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), dimension(:), pointer, contiguous, intent(in) :: ibound
    ! -- local
    ! -- formats
    !
    ! -- adv pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
  end subroutine adv_ar

  !> @brief  Calculate maximum time step length
  !!
  !!  Return the largest time step that meets stability constraints
  !<
  subroutine adv_dt(this, dtmax, msg, thetam)
    ! dummy
    class(TspAdvType) :: this !< this instance
    real(DP), intent(out) :: dtmax !< maximum allowable dt subject to stability constraint
    character(len=*), intent(inout) :: msg !< package/cell dt constraint message
    real(DP), dimension(:), intent(in) :: thetam !< porosity
    ! local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    integer(I4B) :: nrmax
    character(len=LINELENGTH) :: cellstr
    real(DP) :: dt
    real(DP) :: flowmax
    real(DP) :: flowsumpos
    real(DP) :: flowsumneg
    real(DP) :: flownm
    real(DP) :: cell_volume
    dtmax = DNODATA
    nrmax = 0
    msg = ''

    ! If ats_percel not specified by user, then return without making
    ! the courant time step calculation
    if (this%ats_percel == DNODATA) then
      return
    end if

    ! Calculate time step lengths based on stability constraint for each cell
    ! and store the smallest one
    do n = 1, this%dis%nodes
      if (this%ibound(n) == 0) cycle
      flowsumneg = DZERO
      flowsumpos = DZERO
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        flownm = this%fmi%gwfflowja(ipos)
        if (flownm < DZERO) then
          flowsumneg = flowsumneg - flownm
        else
          flowsumpos = flowsumpos + flownm
        end if
      end do
      flowmax = max(flowsumneg, flowsumpos)
      if (flowmax < DPREC) cycle
      cell_volume = this%dis%get_cell_volume(n, this%dis%top(n))
      dt = cell_volume * this%fmi%gwfsat(n) * thetam(n) / flowmax
      dt = dt * this%ats_percel
      if (dt < dtmax) then
        dtmax = dt
        nrmax = n
      end if
    end do
    if (nrmax > 0) then
      call this%dis%noder_to_string(nrmax, cellstr)
      write (msg, *) adjustl(trim(this%memoryPath))//'-'//trim(cellstr)
    end if
  end subroutine adv_dt

  !> @brief  Fill coefficient method for ADV package
  !!
  !!  Method to calculate coefficients and fill amat and rhs.
  !<
  subroutine adv_fc(this, nodes, matrix_sln, idxglo, cnew, rhs)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    integer, intent(in) :: nodes
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    integer(I4B) :: n, m, idiag, ipos
    real(DP) :: omega, qnm
    !
    ! -- Calculate advection terms and add to solution rhs and hcof.  qnm
    !    is the volumetric flow rate and has dimensions of L^/T.
    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      idiag = this%dis%con%ia(n)
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac
        omega = this%adv_weight(this%iadvwt, ipos, n, m, qnm)
        call matrix_sln%add_value_pos(idxglo(ipos), qnm * (DONE - omega))
        call matrix_sln%add_value_pos(idxglo(idiag), qnm * omega)
      end do
    end do
    !
    ! -- TVD
    if (this%iadvwt >= 2) then
      do n = 1, nodes
        if (this%ibound(n) == 0) cycle
        call this%advtvd(n, cnew, rhs)
      end do
    end if
  end subroutine adv_fc

  !> @brief  Calculate TVD
  !!
  !! Use explicit scheme to calculate the advective component of transport.
  !! TVD is an acronym for Total-Variation Diminishing
  !<
  subroutine advtvd(this, n, cnew, rhs)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    ! -- local
    real(DP) :: qtvd
    integer(I4B) :: m, ipos
    !
    ! -- Loop through each n connection.  This will
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      if (this%dis%con%mask(ipos) == 0) cycle
      m = this%dis%con%ja(ipos)
      if (m > n .and. this%ibound(m) /= 0) then
        ! qtvd = this%advqtvd(n, m, ipos, cnew)
        qtvd = this%advqtvd_experimental(n, m, ipos, cnew)
        rhs(n) = rhs(n) - qtvd
        rhs(m) = rhs(m) + qtvd
      end if
    end do
  end subroutine advtvd

  !> @brief  Calculate TVD
  !!
  !! Use explicit scheme to calculate the advective component of transport.
  !! TVD is an acronym for Total-Variation Diminishing
  !<
  function advqtvd(this, n, m, iposnm, cnew) result(qtvd)
    ! -- modules
    use ConstantsModule, only: DPREC
    ! -- return
    real(DP) :: qtvd
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), dimension(:), intent(in) :: cnew
    ! -- local
    integer(I4B) :: ipos, isympos, iup, idn, i2up, j
    real(DP) :: qnm, qmax, qupj, elupdn, elup2up
    real(DP) :: smooth, cdiff, alimiter
    integer(I4B) :: ihc
    integer(I4B) :: ihc2
    real(DP) :: nx, ny, nz
    real(DP) :: nx2, ny2, nz2
    !
    ! -- initialize
    qtvd = DZERO
    !
    ! -- Find upstream node
    isympos = this%dis%con%jas(iposnm)
    qnm = this%fmi%gwfflowja(iposnm)
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n
    else
      iup = n
      idn = m
    end if
    elupdn = this%dis%con%cl1(isympos) + this%dis%con%cl2(isympos)
    !
    ! -- Find connection direction
    ihc = this%dis%con%ihc(isympos)
    call this%dis%connection_normal(n, m, ihc, nx, ny, nz, ipos)
    !
    ! -- Find second node upstream to iup
    i2up = 0
    qmax = DZERO
    do ipos = this%dis%con%ia(iup) + 1, this%dis%con%ia(iup + 1) - 1
      j = this%dis%con%ja(ipos)
      isympos = this%dis%con%jas(ipos)
      ihc2 = this%dis%con%ihc(isympos)
      call this%dis%connection_normal(iup, j, ihc, nx2, ny2, nz2, ipos)

      if (((abs(nx) - abs(nx2)) > 1e-5) .or. (abs(ny) - abs(ny2) > 1e-5)) cycle

      if (this%ibound(j) == 0) cycle
      qupj = this%fmi%gwfflowja(ipos)
     
      if (qupj > qmax) then
        qmax = qupj
        i2up = j
        elup2up = this%dis%con%cl1(isympos) + this%dis%con%cl2(isympos)
      end if
    end do
    !
    ! -- Calculate flux limiting term
    if (i2up > 0) then
      smooth = DZERO
      cdiff = ABS(cnew(iup) - cnew(i2up)) ! adjusted code
      ! cdiff = ABS(cnew(idn) - cnew(iup)) ! original code
      if (cdiff > DPREC) then
        smooth = (cnew(idn) - cnew(iup)) / elupdn * &
          elup2up / (cnew(iup) - cnew(i2up)) ! adjusted code
        ! smooth = (cnew(iup) - cnew(i2up)) / elup2up * &
        !          elupdn / (cnew(idn) - cnew(iup)) ! original code
      end if
      if (smooth > DZERO) then
        select case (this%iadvwt)
        case (2) ! van Leer
          ! alimiter = DTWO * smooth / (DONE + smooth)
          alimiter = min((smooth + dabs(smooth)) / (1.0_dp + dabs(smooth)), 2.0_dp)
        case (3) ! Koren
          alimiter = min(2.0_dp*smooth, 1.0_dp/3.0_dp + 2.0_dp/3.0_dp*smooth, 2.0_dp)
        case (4) ! Superbee
          alimiter = max(min(2.0_dp*smooth, 1.0_dp), min(smooth, 2.0_dp))
        case (5) ! van Albada
          alimiter =  (smooth * smooth + smooth) / (smooth * smooth + 1.0_dp)
        case (6) ! Koren modified
          alimiter = min(4.0_dp*smooth*smooth+smooth,1.0_dp/3.0_dp+2.0_dp/3.0_dp*smooth, 2.0_dp)
        CASE DEFAULT
          alimiter = DZERO
        end select
       
        qtvd = DHALF * alimiter * qnm * (cnew(iup) - cnew(i2up)) ! adjusted code
        ! qtvd = DHALF * alimiter * qnm * (cnew(idn) - cnew(iup)) ! original code
        qtvd = qtvd * this%eqnsclfac
      end if
    end if
  end function advqtvd

  function advqtvd_experimental(this, n, m, iposnm, cnew) result(qtvd)
    ! -- return
    real(DP) :: qtvd
    ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: m
    integer(I4B), intent(in) :: iposnm
    real(DP), dimension(:), intent(in) :: cnew
    ! -- local
    integer(I4B) :: isympos, iup, idn, ihc
    real(DP) :: qnm
    real(DP), dimension(2) :: grad_c
    real(DP) :: smooth, alimiter
    real(DP) :: x_dir, y_dir, z_dir, length
    !
    ! -- initialize
    qtvd = DZERO
    !
    ! -- Find upstream node
    isympos = this%dis%con%jas(iposnm)
    ihc = this%dis%con%ihc(isympos)
    qnm = this%fmi%gwfflowja(iposnm)
    if (qnm > DZERO) then
      ! -- positive flow into n means m is upstream
      iup = m
      idn = n
    else
      iup = n
      idn = m
    end if
    !
    ! Return if straddled cells have same value
    if (abs(cnew(idn) - cnew(iup)) < 1e-5_dp) return
    !
    ! -- Compute cell concentration gradient
    call this%compute_cell_gradient(iup, cnew, grad_c)
    !
    ! -- Compute smoothness factor
    call this%dis%connection_vector(iup, idn, .true., 1.0_dp, 1.0_dp, ihc, x_dir, y_dir, z_dir, length)
    smooth = 2.0_dp * (grad_c(1) * x_dir * length + grad_c(2) * y_dir * length) / (cnew(idn) - cnew(iup)) - 1.0_dp
    !
    ! -- Compute limiter
    ! - TVD
    ! alimiter = max(0.0_dp, min((smooth + dabs(smooth)) / (1.0_dp + dabs(smooth)), 2.0_dp))
    ! - Koren
    ! alimiter = max(0.0_dp, min(2.0_dp*smooth, 1.0_dp/3.0_dp + 2.0_dp/3.0_dp*smooth, 2.0_dp))
    ! - Superbee
    alimiter = max(0.0_dp,max(min(2.0_dp*smooth, 1.0_dp), min(smooth, 2.0_dp)))
    ! -- Compute limited flux
    qtvd = DHALF * alimiter * qnm * (cnew(idn) - cnew(iup)) 
    qtvd = qtvd * this%eqnsclfac


  end function advqtvd_experimental

  subroutine compute_cell_gradient(this, n, cnew, grad_c)
     ! -- dummy
    class(TspAdvType) :: this
    integer(I4B), intent(in) :: n
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(2), intent(out) :: grad_c
    ! -- local
    integer(I4B) :: ipos, isympos, ihc, local_pos
    integer(I4B) :: number_connections
    real(DP), dimension(:, :), allocatable :: d
    real(DP), dimension(:, :), allocatable :: d_trans
    real(DP), dimension(:, :), allocatable :: grad_op
    real(DP), dimension(2, 2) :: g
    real(DP), dimension(2, 2) :: g_inv
    integer(I4B) :: m
    real(DP) :: x_dir, y_dir, z_dir, length

    real(DP), dimension(:), allocatable :: dc

    number_connections = this%dis%con%ia(n + 1) - this%dis%con%ia(n) - 1
    if (number_connections == 1) then
      ! If a cell only has 1 neigbour compute the gradient using finite difference
      ! This case can happen if a triangle element is located in a cornor of a square domain
      ! with two sides being domain boundaries
      ipos = this%dis%con%ia(n) + 1
      m = this%dis%con%ja(ipos)
      isympos = this%dis%con%jas(ipos)
      ihc = this%dis%con%ihc(isympos)
      call this%dis%connection_vector(n, m, .true., 1.0_dp, 1.0_dp, ihc, x_dir, y_dir, z_dir, length)
      
      grad_c(1) = (cnew(m) - cnew(n)) / (x_dir * length)
      grad_c(2) = (cnew(m) - cnew(n)) / (y_dir * length)
      return
    end if

    allocate(d(number_connections, 2))
    allocate(d_trans(2, number_connections))
    allocate(grad_op(2, number_connections))

    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      isympos = this%dis%con%jas(ipos)
      ihc = this%dis%con%ihc(isympos)
     
      call this%dis%connection_vector(n, m, .true., 1.0_dp, 1.0_dp, ihc, x_dir, y_dir, z_dir, length)
      d(local_pos, 1) = x_dir * length
      d(local_pos, 2) = y_dir * length

      d_trans(1, local_pos) = x_dir * length
      d_trans(2, local_pos) = y_dir * length
      local_pos = local_pos + 1
    end do

    g = matmul(d_trans, d)
    g_inv = matinv2(g)
    grad_op = matmul(g_inv, d_trans)

    allocate(dc(number_connections))
    local_pos = 1
    do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
      m = this%dis%con%ja(ipos)
      dc(local_pos) = cnew(m) - cnew(n)
      local_pos = local_pos + 1
    end do

    grad_c = matmul(grad_op, dc)

  end subroutine compute_cell_gradient

  pure function matinv2(A) result(B)
    !! https://fortranwiki.org/fortran/show/Matrix+inversion
    !! Performs a direct calculation of the inverse of a 2×2 matrix.
    real(DP), intent(in) :: A(2,2)   !! Matrix
    real(DP)             :: B(2,2)   !! Inverse matrix
    real(DP)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * A(2,2)
    B(2,1) = -detinv * A(2,1)
    B(1,2) = -detinv * A(1,2)
    B(2,2) = +detinv * A(1,1)
  end function

  !> @brief Calculate advection contribution to flowja
  !<
  subroutine adv_cq(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: nodes
    integer(I4B) :: n, m, idiag, ipos
    real(DP) :: omega, qnm
    !
    ! -- Calculate advection and add to flowja. qnm is the volumetric flow
    !    rate and has dimensions of L^/T.
    nodes = this%dis%nodes
    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      idiag = this%dis%con%ia(n)
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) == 0) cycle
        qnm = this%fmi%gwfflowja(ipos) * this%eqnsclfac
        omega = this%adv_weight(this%iadvwt, ipos, n, m, qnm)
        flowja(ipos) = flowja(ipos) + qnm * omega * cnew(n) + &
                       qnm * (DONE - omega) * cnew(m)
      end do
    end do
    !
    ! -- TVD
    if (this%iadvwt >= 2) call this%advtvd_bd(cnew, flowja)
  end subroutine adv_cq

  !> @brief Add TVD contribution to flowja
  subroutine advtvd_bd(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspAdvType) :: this
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), intent(inout) :: flowja
    ! -- local
    real(DP) :: qtvd, qnm
    integer(I4B) :: nodes, n, m, ipos
    !
    nodes = this%dis%nodes
    do n = 1, nodes
      if (this%ibound(n) == 0) cycle
      do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
        m = this%dis%con%ja(ipos)
        if (this%ibound(m) /= 0) then
          qnm = this%fmi%gwfflowja(ipos)
          ! qtvd = this%advqtvd(n, m, ipos, cnew)
          qtvd = this%advqtvd_experimental(n, m, ipos, cnew)
          flowja(ipos) = flowja(ipos) + qtvd
        end if
      end do
    end do
  end subroutine advtvd_bd

  !> @brief Deallocate memory
  !<
  subroutine adv_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspAdvType) :: this
    !
    ! -- Deallocate arrays if package was active
    if (this%inunit > 0) then
    end if
    !
    ! -- nullify pointers
    this%ibound => null()
    !
    ! -- Scalars
    call mem_deallocate(this%iadvwt)
    call mem_deallocate(this%ats_percel)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine adv_da

  !> @brief Allocate scalars specific to the streamflow energy transport (SFE)
  !! package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspAdvType) :: this
    ! -- local
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%iadvwt, 'IADVWT', this%memoryPath)
    call mem_allocate(this%ats_percel, 'ATS_PERCEL', this%memoryPath)
    !
    ! -- Initialize
    this%iadvwt = 0
    this%ats_percel = DNODATA
    !
    ! -- Advection creates an asymmetric coefficient matrix
    this%iasym = 1
  end subroutine allocate_scalars

  !> @brief Read options
  !!
  !! Read the options block
  !<
  subroutine read_options(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    use SimModule, only: store_error
    ! -- dummy
    class(TspAdvType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg, keyword
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    ! -- formats
    character(len=*), parameter :: fmtiadvwt = &
      &"(4x,'ADVECTION WEIGHTING SCHEME HAS BEEN SET TO: ', a)"
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, blockRequired=.false., &
                              supportOpenClose=.true.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(1x,a)') 'PROCESSING ADVECTION OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('SCHEME')
          call this%parser%GetStringCaps(keyword)
          select case (keyword)
          case ('UPSTREAM')
            this%iadvwt = 0
            write (this%iout, fmtiadvwt) 'UPSTREAM'
          case ('CENTRAL')
            this%iadvwt = 1
            write (this%iout, fmtiadvwt) 'CENTRAL'
          case ('TVD')
            this%iadvwt = 2
            write (this%iout, fmtiadvwt) 'TVD'
          case ('KOREN')
            this%iadvwt = 3
            write (this%iout, fmtiadvwt) 'KOREN'
          case ('SUPERBEE')
            this%iadvwt = 4
            write (this%iout, fmtiadvwt) 'SUPERBEE'
          case ('ALBADA')
            this%iadvwt = 5
            write (this%iout, fmtiadvwt) 'ALBADA'
          case ('KORENMOD')
            this%iadvwt = 6
            write (this%iout, fmtiadvwt) 'KORENMOD'
          case default
            write (errmsg, '(a, a)') &
              'Unknown scheme: ', trim(keyword)
            call store_error(errmsg)
            write (errmsg, '(a, a)') &
              'Scheme must be "UPSTREAM", "CENTRAL" or "TVD"'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end select
        case ('ATS_PERCEL')
          this%ats_percel = this%parser%GetDouble()
          if (this%ats_percel == DZERO) this%ats_percel = DNODATA
          write (this%iout, '(4x,a,1pg15.6)') &
            'User-specified fractional cell distance for adaptive time &
            &steps: ', this%ats_percel
        case default
          write (errmsg, '(a,a)') 'Unknown ADVECTION option: ', &
            trim(keyword)
          call store_error(errmsg, terminate=.TRUE.)
        end select
      end do
      write (this%iout, '(1x,a)') 'END OF ADVECTION OPTIONS'
    end if
  end subroutine read_options

  !> @ brief Advection weight
  !!
  !! Calculate the advection weight
  !<
  function adv_weight(this, iadvwt, ipos, n, m, qnm) result(omega)
    ! -- return
    real(DP) :: omega
    ! -- dummy
    class(TspAdvType) :: this
    integer, intent(in) :: iadvwt
    integer, intent(in) :: ipos
    integer, intent(in) :: n
    integer, intent(in) :: m
    real(DP), intent(in) :: qnm
    ! -- local
    real(DP) :: lnm, lmn

    select case (iadvwt)
    case (1)
      ! -- calculate weight based on distances between nodes and the shared
      !    face of the connection
      if (this%dis%con%ihc(this%dis%con%jas(ipos)) == 0) then
        ! -- vertical connection; assume cell is fully saturated
        lnm = DHALF * (this%dis%top(n) - this%dis%bot(n))
        lmn = DHALF * (this%dis%top(m) - this%dis%bot(m))
      else
        ! -- horizontal connection
        lnm = this%dis%con%cl1(this%dis%con%jas(ipos))
        lmn = this%dis%con%cl2(this%dis%con%jas(ipos))
      end if
      omega = lmn / (lnm + lmn)
    case (0, 2 :)
      ! -- use upstream weighting for upstream and tvd schemes
      if (qnm > DZERO) then
        omega = DZERO
      else
        omega = DONE
      end if
    end select
  end function adv_weight

end module TspAdvModule
