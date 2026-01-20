!> @brief This module contains the storage package methods
!!
!! This module contains the methods used to add the effects of storage
!! on the surface water flow equation.
!!
!<
module SwfStoModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, LENBUDTXT, LINELENGTH, LENMEMPATH
  use MemoryHelperModule, only: create_mem_path
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, store_error_filename
  use BaseDisModule, only: DisBaseType
  use NumericalPackageModule, only: NumericalPackageType
  use MatrixBaseModule
  use Disv1dModule, only: Disv1dType
  use SwfCxsModule, only: SwfCxsType

  implicit none
  public :: SwfStoType, sto_cr

  character(len=LENBUDTXT), dimension(1) :: budtxt = & !< text labels for budget terms
    &['         STORAGE']

  type, extends(NumericalPackageType) :: SwfStoType
    integer(I4B), pointer :: iss => null() !< steady state flag: 1 = steady, 0 = transient
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    real(DP), dimension(:), pointer, contiguous :: qsto => null() !< storage rates

    ! -- pointers to information in dfw package
    integer(I4B), dimension(:), pointer, contiguous :: idcxs => null() !< pointer to cross section id vector in dfw

    ! -- pointer to packages needed for calculations
    type(SwfCxsType), pointer :: cxs

    ! -- input context pointers for read and prepare
    integer(I4B), pointer :: iper => null() !< input context loaded period
    character(len=:), pointer :: storage !< input context storage string
  contains
    procedure :: sto_ar
    procedure :: sto_rp
    procedure :: sto_ad
    procedure :: sto_fc
    procedure :: sto_fc_dis1d
    procedure :: sto_fc_dis2d
    !procedure :: sto_fn
    procedure :: sto_cq
    procedure :: sto_bd
    procedure :: sto_save_model_flows
    procedure :: sto_da
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_data
    procedure, private :: log_options
    procedure, private :: set_dfw_pointers
    procedure, private :: reach_length_pointer
    procedure, private :: calc_storage_dis1d
    procedure, private :: calc_storage_dis2d
  end type

contains

  !> @ brief Create a new package object
  !!
  !!  Create a new storage (STO) object
  !!
  !<
  subroutine sto_cr(stoobj, name_model, mempath, inunit, iout, cxs)
    ! -- dummy variables
    type(SwfStoType), pointer :: stoobj !< SwfStoType object
    character(len=*), intent(in) :: name_model !< name of model
    character(len=*), intent(in) :: mempath !< input context mem path
    integer(I4B), intent(in) :: inunit !< package input file unit
    integer(I4B), intent(in) :: iout !< model listing file unit
    type(SwfCxsType), pointer, intent(in) :: cxs !< the pointer to the cxs package
    !
    ! -- Create the object
    allocate (stoobj)
    !
    ! -- create name and memory path
    call stoobj%set_names(1, name_model, 'STO', 'STO', mempath)
    !
    ! -- Allocate scalars
    call stoobj%allocate_scalars()
    !
    ! -- Set variables
    stoobj%inunit = inunit
    stoobj%iout = iout

    ! -- store pointers
    stoobj%cxs => cxs

  end subroutine sto_cr

  !> @ brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the STO package.
  !!
  !<
  subroutine sto_ar(this, dis, ibound)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization object
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< model ibound array
    ! -- local variables
    ! -- formats
    character(len=*), parameter :: fmtsto = &
      "(1x,/1x,'STO -- STORAGE PACKAGE, VERSION 1, 10/27/2023', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! --print a message identifying the storage package.
    write (this%iout, fmtsto) this%inunit

    ! -- set pointers to data in dfw package
    call this%set_dfw_pointers()

    !
    ! -- store pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
    !
    ! -- set pointer to model iss
    call mem_setptr(this%iss, 'ISS', create_mem_path(this%name_model))
    !
    ! -- Allocate arrays
    call this%allocate_arrays(dis%nodes)
    !
    ! -- Read storage options
    call this%source_options()
    !
    ! -- read the data block
    ! no griddata at the moment for SWF Storage Package
    ! call this%source_data()
  end subroutine sto_ar

  !> @ brief Read and prepare method for package
  !!
  !!  Method to read and prepare stress period data for the STO package.
  !!
  !<
  subroutine sto_rp(this)
    ! -- modules
    use TdisModule, only: kper
    implicit none
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    ! -- local variables
    character(len=16) :: css(0:1)
    ! -- data
    data css(0)/'       TRANSIENT'/
    data css(1)/'    STEADY-STATE'/
    !
    ! -- confirm package is active
    if (this%inunit <= 0) return
    !
    ! -- confirm loaded iper
    if (this%iper /= kper) return
    !
    write (this%iout, '(//,1x,a)') 'PROCESSING STORAGE PERIOD DATA'
    !
    ! -- set period iss
    if (this%storage == 'STEADY-STATE') then
      this%iss = 1
    else if (this%storage == 'TRANSIENT') then
      this%iss = 0
    else
      write (errmsg, '(a,a)') 'Unknown STORAGE data tag: ', &
        trim(this%storage)
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
    !
    write (this%iout, '(1x,a)') 'END PROCESSING STORAGE PERIOD DATA'
    !
    write (this%iout, '(//1X,A,I0,A,A,/)') &
      'STRESS PERIOD ', kper, ' IS ', trim(adjustl(css(this%iss)))
  end subroutine sto_rp

  !> @ brief Advance the package
  !!
  !!  Advance data in the STO package.
  !!
  !<
  subroutine sto_ad(this)
    ! -- modules
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
  end subroutine sto_ad

  !> @ brief Fill A and right-hand side for the package
  !!
  !!  Fill the coefficient matrix and right-hand side with the STO package terms.
  !!
  !<
  subroutine sto_fc(this, kiter, stage_old, stage_new, matrix_sln, idxglo, rhs)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(SwfStoType) :: this
    integer(I4B) :: kiter
    real(DP), intent(inout), dimension(:) :: stage_old
    real(DP), intent(inout), dimension(:) :: stage_new
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    ! -- local
    ! -- formats
    character(len=*), parameter :: fmtsperror = &
      &"('Detected time step length of zero.  SWF Storage Package cannot be ', &
      &'used unless delt is non-zero.')"
    !
    ! -- test if steady-state stress period
    if (this%iss /= 0) return
    !
    ! -- Ensure time step length is not zero
    if (delt == DZERO) then
      write (errmsg, fmtsperror)
      call store_error(errmsg, terminate=.TRUE.)
    end if

    if (this%dis%is_1d()) then
      call this%sto_fc_dis1d(kiter, stage_old, stage_new, matrix_sln, idxglo, rhs)
    else
      call this%sto_fc_dis2d(kiter, stage_old, stage_new, matrix_sln, idxglo, rhs)
    end if

  end subroutine sto_fc

  !> @ brief Fill A and right-hand side for the package
  !!
  !!  Fill the coefficient matrix and right-hand side with the STO package terms.
  !!
  !<
  subroutine sto_fc_dis1d(this, kiter, stage_old, stage_new, matrix_sln, &
                          idxglo, rhs)
    ! -- modules
    ! -- dummy
    class(SwfStoType) :: this
    integer(I4B) :: kiter
    real(DP), intent(inout), dimension(:) :: stage_old
    real(DP), intent(inout), dimension(:) :: stage_new
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: derv
    real(DP) :: qsto
    real(DP), dimension(:), pointer :: reach_length

    ! Set pointer to reach_length for 1d
    reach_length => this%reach_length_pointer()

    ! -- Calculate coefficients and put into amat
    do n = 1, this%dis%nodes

      ! -- skip if constant stage
      if (this%ibound(n) < 0) cycle

      call this%calc_storage_dis1d(n, stage_new(n), stage_old(n), &
                                   reach_length(n), qsto, derv)

      ! -- Fill amat and rhs
      idiag = this%dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), -derv)
      rhs(n) = rhs(n) + qsto - derv * stage_new(n)

    end do
  end subroutine sto_fc_dis1d

  !> @ brief Fill A and right-hand side for the package
  !!
  !!  Fill the coefficient matrix and right-hand side with the STO package terms.
  !!
  !<
  subroutine sto_fc_dis2d(this, kiter, stage_old, stage_new, matrix_sln, &
                          idxglo, rhs)
    ! -- modules
    ! -- dummy
    class(SwfStoType) :: this
    integer(I4B) :: kiter
    real(DP), intent(inout), dimension(:) :: stage_old
    real(DP), intent(inout), dimension(:) :: stage_new
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    ! -- local
    integer(I4B) :: n, idiag
    real(DP) :: derv
    real(DP) :: qsto

    ! -- Calculate coefficients and put into amat
    do n = 1, this%dis%nodes
      !
      ! -- skip if constant stage
      if (this%ibound(n) < 0) cycle

      ! Calculate storage and derivative term
      call this%calc_storage_dis2d(n, stage_new(n), stage_old(n), &
                                   qsto, derv)

      ! -- Fill amat and rhs
      idiag = this%dis%con%ia(n)
      call matrix_sln%add_value_pos(idxglo(idiag), -derv)
      rhs(n) = rhs(n) + qsto - derv * stage_new(n)

    end do

  end subroutine sto_fc_dis2d

  !> @ brief Calculate flows for package
  !<
  subroutine sto_cq(this, flowja, stage_new, stage_old)
    ! -- dummy
    class(SwfStoType) :: this
    real(DP), intent(inout), dimension(:) :: flowja
    real(DP), intent(inout), dimension(:) :: stage_new
    real(DP), intent(inout), dimension(:) :: stage_old
    ! -- local
    real(DP), dimension(:), pointer :: reach_length
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: dx
    real(DP) :: q

    ! -- test if steady-state stress period
    if (this%iss /= 0) return

    ! Set pointer to reach_length for 1d
    reach_length => this%reach_length_pointer()

    ! -- Calculate storage term
    do n = 1, this%dis%nodes
      !
      ! -- skip if constant stage
      if (this%ibound(n) < 0) cycle

      ! Calculate storage for either the DIS1D or DIS2D cases and
      ! add to flowja
      if (associated(reach_length)) then
        dx = reach_length(n)
        call this%calc_storage_dis1d(n, stage_new(n), stage_old(n), dx, q)
      else
        call this%calc_storage_dis2d(n, stage_new(n), stage_old(n), q)
      end if
      this%qsto(n) = -q
      idiag = this%dis%con%ia(n)
      flowja(idiag) = flowja(idiag) + this%qsto(n)

    end do
  end subroutine sto_cq

  subroutine calc_storage_dis1d(this, n, stage_new, stage_old, dx, qsto, derv)
    ! module
    use TdisModule, only: delt
    use MathUtilModule, only: get_perturbation
    ! dummy
    class(SwfStoType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: stage_new
    real(DP), intent(in) :: stage_old
    real(DP), intent(in) :: dx
    real(DP), intent(inout) :: qsto
    real(DP), intent(inout), optional :: derv
    ! local
    real(DP) :: depth_new
    real(DP) :: depth_old
    real(DP) :: width_n
    real(DP) :: width_m
    real(DP) :: cxs_area_new
    real(DP) :: cxs_area_old
    real(DP) :: cxs_area_eps
    real(DP) :: eps

    call this%dis%get_flow_width(n, n, 0, width_n, width_m)
    depth_new = stage_new - this%dis%bot(n)
    depth_old = stage_old - this%dis%bot(n)
    cxs_area_new = this%cxs%get_area(this%idcxs(n), width_n, depth_new)
    cxs_area_old = this%cxs%get_area(this%idcxs(n), width_n, depth_old)
    qsto = (cxs_area_new - cxs_area_old) * dx / delt
    if (present(derv)) then
      eps = get_perturbation(depth_new)
      cxs_area_eps = this%cxs%get_area(this%idcxs(n), width_n, depth_new + eps)
      derv = (cxs_area_eps - cxs_area_new) * dx / delt / eps
    end if

  end subroutine calc_storage_dis1d

  subroutine calc_storage_dis2d(this, n, stage_new, stage_old, qsto, derv)
    ! module
    use TdisModule, only: delt
    use MathUtilModule, only: get_perturbation
    ! dummy
    class(SwfStoType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: stage_new
    real(DP), intent(in) :: stage_old
    real(DP), intent(inout) :: qsto
    real(DP), intent(inout), optional :: derv
    ! local
    real(DP) :: area
    real(DP) :: depth_new
    real(DP) :: depth_old
    real(DP) :: depth_eps
    real(DP) :: volume_new
    real(DP) :: volume_old
    real(DP) :: eps

    area = this%dis%get_area(n)
    depth_new = stage_new - this%dis%bot(n)
    depth_old = stage_old - this%dis%bot(n)
    volume_new = area * depth_new
    volume_old = area * depth_old
    qsto = (volume_new - volume_old) / delt

    if (present(derv)) then
      eps = get_perturbation(depth_new)
      depth_eps = depth_new + eps
      derv = (depth_eps - depth_new) * area / delt / eps
    end if

  end subroutine calc_storage_dis2d

  !> @ brief Model budget calculation for package
  !!
  !!  Budget calculation for the STO package components. Components include
  !!  specific storage and specific yield storage.
  !!
  !<
  subroutine sto_bd(this, isuppress_output, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    integer(I4B), intent(in) :: isuppress_output !< flag to suppress model output
    type(BudgetType), intent(inout) :: model_budget !< model budget object
    ! -- local variables
    real(DP) :: rin
    real(DP) :: rout
    !
    ! -- Add storage rates to model budget
    call rate_accumulator(this%qsto, rin, rout)
    call model_budget%addentry(rin, rout, delt, '             STO', &
                               isuppress_output, '         STORAGE')
  end subroutine sto_bd

  !> @ brief Save model flows for package
  !!
  !!  Save cell-by-cell budget terms for the STO package.
  !!
  !<
  subroutine sto_save_model_flows(this, icbcfl, icbcun)
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    integer(I4B), intent(in) :: icbcfl !< flag to output budget data
    integer(I4B), intent(in) :: icbcun !< cell-by-cell file unit number
    ! -- local variables
    integer(I4B) :: ibinun
    integer(I4B) :: iprint, nvaluesp, nwidthp
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    real(DP) :: dinact
    !
    ! -- Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! -- Record the storage rates if requested
    if (ibinun /= 0) then
      iprint = 0
      dinact = DZERO
      !
      ! -- qsto
      call this%dis%record_array(this%qsto, this%iout, iprint, -ibinun, &
                                 budtxt(1), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)
    end if
  end subroutine sto_save_model_flows

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate STO package scalars and arrays.
  !!
  !<
  subroutine sto_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    !
    ! -- Deallocate arrays if package is active
    if (this%inunit > 0) then
      call mem_deallocate(this%qsto)
      nullify (this%idcxs)
      nullify (this%iper)
      nullify (this%storage)
    end if
    !
    ! -- Deallocate scalars
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine sto_da

  !> @ brief Allocate scalars
  !!
  !! Allocate and initialize scalars for the STO package. The base numerical
  !! package allocate scalars method is also called.
  !!
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- allocate scalars
    !call mem_allocate(this%xxx, 'XXX', this%memoryPath)
    !
    ! -- initialize scalars
    !this%xxx = 0
  end subroutine allocate_scalars

  !> @ brief Allocate package arrays
  !!
  !!  Allocate and initialize STO package arrays.
  !!
  !<
  subroutine allocate_arrays(this, nodes)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy variables
    class(SwfStoType), target :: this !< SwfStoType object
    integer(I4B), intent(in) :: nodes !< active model nodes
    ! -- local variables
    integer(I4B) :: n
    !
    ! -- Allocate arrays
    call mem_allocate(this%qsto, nodes, 'STRGSS', this%memoryPath)
    !
    ! -- set input context pointers
    if (this%inunit > 0) then
      call mem_setptr(this%iper, 'IPER', this%input_mempath)
      call mem_setptr(this%storage, 'STORAGE', this%input_mempath)
    end if
    !
    ! -- Initialize arrays
    this%iss = 0
    do n = 1, nodes
      this%qsto(n) = DZERO
    end do
  end subroutine allocate_arrays

  !> @ brief Source input options for package
  !!
  !!  Source options block parameters for STO package.
  !!
  !<
  subroutine source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use SourceCommonModule, only: filein_fname
    use SwfStoInputModule, only: SwfStoParamFoundType
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    ! -- local variables
    type(SwfStoParamFoundType) :: found
    !
    ! -- source package input
    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)
    !
    if (found%ipakcb) then
      this%ipakcb = -1
    end if
    !
    ! -- log found options
    call this%log_options(found)
  end subroutine source_options

  !> @ brief Log found options for package
  !!
  !!  Log options block for STO package.
  !!
  !<
  subroutine log_options(this, found)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use SwfStoInputModule, only: SwfStoParamFoundType
    ! -- dummy variables
    class(SwfStoType) :: this !< SwfStoType object
    type(SwfStoParamFoundType), intent(in) :: found
    ! -- local variables
    ! -- formats
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO.')"
    !
    write (this%iout, '(1x,a)') 'PROCESSING STORAGE OPTIONS'
    !
    if (found%ipakcb) then
      write (this%iout, fmtisvflow)
    end if
    !
    write (this%iout, '(1x,a)') 'END OF STORAGE OPTIONS'
  end subroutine log_options

  !> @ brief Source input data for package
  !!
  !!  Source griddata block parameters for STO package.
  !!
  !<
  subroutine source_data(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use SwfStoInputModule, only: SwfStoParamFoundType
    ! -- dummy variables
    class(SwfStotype) :: this !< SwfStoType object
    ! -- local variables
    character(len=24), dimension(1) :: aname
    integer(I4B), dimension(:), pointer, contiguous :: map
    !type(SwfStoParamFoundType) :: found
    !
    ! -- initialize data
    data aname(1)/'                     XXX'/
    !
    ! -- set map to reduce data input arrays
    map => null()
    if (this%dis%nodes < this%dis%nodesuser) map => this%dis%nodeuser
    !
    ! -- log griddata
    write (this%iout, '(1x,a)') 'PROCESSING GRIDDATA'
    write (this%iout, '(1x,a)') 'END PROCESSING GRIDDATA'
  end subroutine source_data

  !> @brief Set pointers to channel properties in DFW Package
  !<
  subroutine set_dfw_pointers(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr, mem_allocate
    ! -- dummy
    class(SwfStoType) :: this !< this instance
    ! -- local
    character(len=LENMEMPATH) :: dfw_mem_path

    dfw_mem_path = create_mem_path(this%name_model, 'DFW')
    call mem_setptr(this%idcxs, 'IDCXS', dfw_mem_path)

  end subroutine set_dfw_pointers

  function reach_length_pointer(this) result(ptr)
    ! dummy
    class(SwfStoType) :: this !< this instance
    ! return
    real(DP), dimension(:), pointer :: ptr
    ! local
    class(DisBaseType), pointer :: dis

    ptr => null()
    dis => this%dis
    select type (dis)
    type is (Disv1dType)
      ptr => dis%length
    end select

  end function reach_length_pointer

end module SwfStoModule
