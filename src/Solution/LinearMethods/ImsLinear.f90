module IMSLinearModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENSOLUTIONNAME, LENMEMPATH, &
                             IZERO, DZERO, DPREC, DSAME, &
                             DEM8, DEM6, DEM5, DEM4, DEM3, DEM2, DEM1, &
                             DHALF, DONE, DTWO, &
                             VDEBUG
  use IMSLinearBaseModule, only: ims_base_cg, ims_base_bcgs, &
                                 ims_base_calc_order, &
                                 ims_base_scale, ims_base_pcu, &
                                 ims_base_residual, ims_base_epfact, &
                                 ims_calc_pcdims
  use IMSLinearMisc, only: ims_base_pccrs
  use BlockParserModule, only: BlockParserType
  use MatrixBaseModule
  use ConvergenceSummaryModule
  use ImsLinearSettingsModule
  use ImsLinearAmgModule, only: ImsAmgDataType, ims_amg_setup, ims_amg_da, &
                                ims_amg_summary
  use ProfilerModule, only: g_prof

  implicit none
  private

  type, public :: ImsLinearDataType
    character(len=LENMEMPATH) :: memoryPath !< the path for storing variables in the memory manager
    integer(I4B), pointer :: iout => null() !< simulation listing file unit
    integer(I4B), pointer :: IPRIMS => null() !< print flag
    ! input variables (pointing to fields in input structure)
    real(DP), pointer :: DVCLOSE => null() !< dependent variable closure criterion
    real(DP), pointer :: RCLOSE => null() !< residual closure criterion
    integer(I4B), pointer :: ICNVGOPT => null() !< convergence option
    integer(I4B), pointer :: ITER1 => null() !< max. iterations
    integer(I4B), pointer :: ILINMETH => null() !< linear solver method
    integer(I4B), pointer :: ISCL => null() !< scaling method
    integer(I4B), pointer :: IORD => null() !< reordering method
    integer(I4B), pointer :: NORTH => null() !< number of orthogonalizations

    real(DP), pointer :: RELAX => null() !< relaxation factor
    integer(I4B), pointer :: LEVEL => null() !< nr. of preconditioner levels
    real(DP), pointer :: DROPTOL => null() !< drop tolerance for preconditioner
    integer(I4B), pointer :: NSMOOTH => null() !< AMG smoother iterations per level
    integer(I4B), pointer :: ISMOOTHER => null() !< AMG smoother type
    real(DP), pointer :: STHRESH => null() !< AMG strength-of-connection threshold
    !
    integer(I4B), pointer :: IPC => null() !< preconditioner flag
    integer(I4B), pointer :: IACPC => null() !< preconditioner CRS row pointers
    integer(I4B), pointer :: NITERC => null() !<
    integer(I4B), pointer :: NIABCGS => null() !< size of working vectors for BCGS linear accelerator
    integer(I4B), pointer :: NIAPC => null() !< preconditioner number of rows
    integer(I4B), pointer :: NJAPC => null() !< preconditioner number of non-zero entries
    real(DP), pointer :: EPFACT => null() !< factor for decreasing convergence criteria in seubsequent Picard iterations
    real(DP), pointer :: L2NORM0 => null() !< initial L2 norm
    ! -- ilut variables
    integer(I4B), pointer :: NJLU => null() !< length of jlu work vector
    integer(I4B), pointer :: NJW => null() !< length of jw work vector
    integer(I4B), pointer :: NWLU => null() !< length of wlu work vector
    ! -- pointers to solution variables
    integer(I4B), pointer :: NEQ => null() !< number of equations (rows in matrix)
    integer(I4B), pointer :: NJA => null() !< number of non-zero values in amat
    integer(I4B), dimension(:), pointer, contiguous :: IA => null() !< position of start of each row
    integer(I4B), dimension(:), pointer, contiguous :: JA => null() !< column pointer
    real(DP), dimension(:), pointer, contiguous :: AMAT => null() !< coefficient matrix
    real(DP), dimension(:), pointer, contiguous :: RHS => null() !< right-hand side of equation
    real(DP), dimension(:), pointer, contiguous :: X => null() !< dependent variable
    ! VECTORS
    real(DP), pointer, dimension(:), contiguous :: DSCALE => null() !< scaling factor
    real(DP), pointer, dimension(:), contiguous :: DSCALE2 => null() !< unscaling factor
    integer(I4B), pointer, dimension(:), contiguous :: IAPC => null() !< position of start of each row in preconditioner matrix
    integer(I4B), pointer, dimension(:), contiguous :: JAPC => null() !< preconditioner matrix column pointer
    real(DP), pointer, dimension(:), contiguous :: APC => null() !< preconditioner coefficient matrix
    integer(I4B), pointer, dimension(:), contiguous :: LORDER => null() !< reordering mapping
    integer(I4B), pointer, dimension(:), contiguous :: IORDER => null() !< mapping to restore reordered matrix
    integer(I4B), pointer, dimension(:), contiguous :: IARO => null() !< position of start of each row in reordered matrix
    integer(I4B), pointer, dimension(:), contiguous :: JARO => null() !< reordered matrix column pointer
    real(DP), pointer, dimension(:), contiguous :: ARO => null() !< reordered coefficient matrix
    ! WORKING ARRAYS
    integer(I4B), pointer, dimension(:), contiguous :: IW => null() !< integer working array
    real(DP), pointer, dimension(:), contiguous :: W => null() !< real working array
    integer(I4B), pointer, dimension(:), contiguous :: ID => null() !< integer working array
    real(DP), pointer, dimension(:), contiguous :: D => null() !< real working array
    real(DP), pointer, dimension(:), contiguous :: P => null() !< real working array
    real(DP), pointer, dimension(:), contiguous :: Q => null() !< real working array
    real(DP), pointer, dimension(:), contiguous :: Z => null() !< real working array
    ! BICGSTAB WORKING ARRAYS
    real(DP), pointer, dimension(:), contiguous :: T => null() !< BICGSTAB real working array
    real(DP), pointer, dimension(:), contiguous :: V => null() !< BICGSTAB real working array
    real(DP), pointer, dimension(:), contiguous :: DHAT => null() !< BICGSTAB real working array
    real(DP), pointer, dimension(:), contiguous :: PHAT => null() !< BICGSTAB real working array
    real(DP), pointer, dimension(:), contiguous :: QHAT => null() !< rBICGSTAB eal working array
    ! POINTERS FOR USE WITH BOTH ORIGINAL AND RCM ORDERINGS
    integer(I4B), pointer, dimension(:), contiguous :: IA0 => null() !< pointer to current CRS row pointers
    integer(I4B), pointer, dimension(:), contiguous :: JA0 => null() !< pointer to current CRS column pointers
    real(DP), pointer, dimension(:), contiguous :: A0 => null() !< pointer to current coefficient matrix
    ! ILUT WORKING ARRAYS
    integer(I4B), pointer, dimension(:), contiguous :: JLU => null() !< ilut integer working array
    integer(I4B), pointer, dimension(:), contiguous :: JW => null() !< ilut integer working array
    real(DP), pointer, dimension(:), contiguous :: WLU => null() !< ilut real working array

    ! AMG preconditioner (persistent across outer iterations)
    type(ImsAmgDataType) :: amg_prec
    ! PROCEDURES (METHODS)
  contains
    procedure :: IMSLINEAR_ALLOCATE => imslinear_ar
    procedure :: imslinear_summary
    procedure :: IMSLINEAR_APPLY => imslinear_ap
    procedure :: IMSLINEAR_DA => imslinear_da
    procedure, private :: allocate_scalars
    ! -- PRIVATE PROCEDURES
    procedure, private :: SET_IMSLINEAR_INPUT => imslinear_set_input
  end type ImsLinearDataType

contains

  !> @ brief Allocate storage and read data
    !!
    !!  Allocate storage for linear accelerators and read data
    !!
  !<
  subroutine imslinear_ar(this, NAME, IOUT, IPRIMS, MXITER, &
                          NEQ, matrix, RHS, X, linear_settings)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use MemoryHelperModule, only: create_mem_path
    use SimModule, only: store_error, count_errors, &
                         deprecation_warning
    ! -- dummy variables
    class(ImsLinearDataType), intent(INOUT) :: this !< ImsLinearDataType instance
    character(LEN=LENSOLUTIONNAME), intent(IN) :: NAME !< solution name
    integer(I4B), intent(IN) :: IOUT !< simulation listing file unit
    integer(I4B), target, intent(IN) :: IPRIMS !< print option
    integer(I4B), intent(IN) :: MXITER !< maximum outer iterations
    integer(I4B), target, intent(IN) :: NEQ !< number of equations
    class(MatrixBaseType), pointer :: matrix
    real(DP), dimension(NEQ), target, intent(INOUT) :: RHS !< right-hand side
    real(DP), dimension(NEQ), target, intent(INOUT) :: X !< dependent variables
    type(ImsLinearSettingsType), pointer :: linear_settings !< the settings form the IMS file
    ! -- local variables
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: n
    integer(I4B) :: i0
    integer(I4B) :: iscllen, iolen

    !
    ! -- DEFINE NAME
    this%memoryPath = create_mem_path(name, 'IMSLINEAR')
    !
    ! -- SET pointers to IMS settings
    this%DVCLOSE => linear_settings%dvclose
    this%RCLOSE => linear_settings%rclose
    this%ICNVGOPT => linear_settings%icnvgopt
    this%ITER1 => linear_settings%iter1
    this%ILINMETH => linear_settings%ilinmeth
    this%ISCL => linear_settings%iscl
    this%IORD => linear_settings%iord
    this%NORTH => linear_settings%north
    this%RELAX => linear_settings%relax
    this%LEVEL => linear_settings%level
    this%DROPTOL => linear_settings%droptol
    this%NSMOOTH => linear_settings%nsmooth
    this%ISMOOTHER => linear_settings%smoother_type
    this%STHRESH => linear_settings%strength_threshold
    !
    ! -- SET POINTERS TO SOLUTION STORAGE
    this%IPRIMS => IPRIMS
    this%NEQ => NEQ
    call matrix%get_aij(this%IA, this%JA, this%AMAT)
    call mem_allocate(this%NJA, 'NJA', this%memoryPath)
    this%NJA = size(this%AMAT)
    this%RHS => RHS
    this%X => X
    !
    ! -- ALLOCATE SCALAR VARIABLES
    call this%allocate_scalars()
    !
    ! -- initialize iout
    this%iout = iout
    !
    ! -- DEFAULT VALUES
    this%IPC = IPC_UNKNOWN
    !
    this%IACPC = 0
    !
    ! -- PRINT A MESSAGE IDENTIFYING IMSLINEAR SOLVER PACKAGE
    write (iout, 2000)
02000 format(1x, /1x, 'IMSLINEAR -- UNSTRUCTURED LINEAR SOLUTION', &
           ' PACKAGE, VERSION 8, 04/28/2017')
    !
    ! -- DETERMINE PRECONDITIONER
    this%IPC = resolve_ipc(linear_settings%ipc_type, this%LEVEL, this%RELAX)
    !
    ! -- ERROR CHECKING FOR OPTIONS
    if (this%ISCL < 0) this%ISCL = 0
    if (this%ISCL > 2) then
      write (errmsg, '(A)') 'IMSLINEAR7AR ISCL MUST BE <= 2'
      call store_error(errmsg)
    end if
    if (this%IORD < 0) this%IORD = 0
    if (this%IORD > 2) then
      write (errmsg, '(A)') 'IMSLINEAR7AR IORD MUST BE <= 2'
      call store_error(errmsg)
    end if
    if (this%NORTH < 0) then
      write (errmsg, '(A)') 'IMSLINEAR7AR NORTH MUST >= 0'
      call store_error(errmsg)
    end if
    if (this%RCLOSE == DZERO) then
      if (this%ICNVGOPT /= CNVG_REL_L2_NORM) then
        write (errmsg, '(A)') 'IMSLINEAR7AR RCLOSE MUST > 0.0'
        call store_error(errmsg)
      end if
    end if
    if (this%RELAX < DZERO) then
      write (errmsg, '(A)') 'IMSLINEAR7AR RELAX MUST BE >= 0.0'
      call store_error(errmsg)
    end if
    if (this%RELAX > DONE) then
      write (errmsg, '(A)') 'IMSLINEAR7AR RELAX MUST BE <= 1.0'
      call store_error(errmsg)
    end if
    if (this%IPC == IPC_AMG .and. this%ILINMETH == BCGS_METHOD .and. &
        this%RELAX /= DONE) then
      write (errmsg, '(A)') 'PRECONDITIONER_TYPE AMG with '// &
        'LINEAR_ACCELERATION BICGSTAB requires RELAXATION_FACTOR = 1. '// &
        'Adaptive omega (RELAXATION_FACTOR 0) makes the preconditioner '// &
        'variable and disrupts BiCGSTAB Krylov optimality; use '// &
        'RELAXATION_FACTOR 1 (fixed undamped) with BICGSTAB and AMG.'
      call store_error(errmsg)
    end if
    !
    ! -- INITIALIZE IMSLINEAR VARIABLES
    this%NITERC = 0
    !
    ! -- ALLOCATE AND INITIALIZE MEMORY FOR IMSLINEAR
    iscllen = 1
    if (this%ISCL /= SCL_NONE) iscllen = NEQ
    call mem_allocate(this%DSCALE, iscllen, 'DSCALE', trim(this%memoryPath))
    call mem_allocate(this%DSCALE2, iscllen, 'DSCALE2', trim(this%memoryPath))
    !
    ! -- allocate and initialize the preconditioner work arrays
    call precond_allocate(this)
    !
    ! -- ALLOCATE SPACE FOR PERMUTATION VECTOR
    i0 = 1
    iolen = 1
    if (this%IORD /= ORD_NONE) then
      i0 = this%NEQ
      iolen = this%NJA
    end if
    call mem_allocate(this%LORDER, i0, 'LORDER', trim(this%memoryPath))
    call mem_allocate(this%IORDER, i0, 'IORDER', trim(this%memoryPath))
    call mem_allocate(this%IARO, i0 + 1, 'IARO', trim(this%memoryPath))
    call mem_allocate(this%JARO, iolen, 'JARO', trim(this%memoryPath))
    call mem_allocate(this%ARO, iolen, 'ARO', trim(this%memoryPath))
    !
    ! -- ALLOCATE WORKING VECTORS FOR IMSLINEAR SOLVER
    call mem_allocate(this%ID, this%NEQ, 'ID', trim(this%memoryPath))
    call mem_allocate(this%D, this%NEQ, 'D', trim(this%memoryPath))
    call mem_allocate(this%P, this%NEQ, 'P', trim(this%memoryPath))
    call mem_allocate(this%Q, this%NEQ, 'Q', trim(this%memoryPath))
    call mem_allocate(this%Z, this%NEQ, 'Z', trim(this%memoryPath))
    !
    ! -- ALLOCATE MEMORY FOR BCGS WORKING ARRAYS
    this%NIABCGS = 1
    if (this%ILINMETH == BCGS_METHOD) then
      this%NIABCGS = this%NEQ
    end if
    call mem_allocate(this%T, this%NIABCGS, 'T', trim(this%memoryPath))
    call mem_allocate(this%V, this%NIABCGS, 'V', trim(this%memoryPath))
    call mem_allocate(this%DHAT, this%NIABCGS, 'DHAT', trim(this%memoryPath))
    call mem_allocate(this%PHAT, this%NIABCGS, 'PHAT', trim(this%memoryPath))
    call mem_allocate(this%QHAT, this%NIABCGS, 'QHAT', trim(this%memoryPath))
    !
    ! -- INITIALIZE IMSLINEAR VECTORS
    do n = 1, iscllen
      this%DSCALE(n) = DONE
      this%DSCALE2(n) = DONE
    end do
    !
    ! -- WORKING VECTORS
    do n = 1, this%NEQ
      this%ID(n) = IZERO
      this%D(n) = DZERO
      this%P(n) = DZERO
      this%Q(n) = DZERO
      this%Z(n) = DZERO
    end do
    do n = 1, this%NIAPC
      this%IW(n) = IZERO
      this%W(n) = DZERO
    end do
    !
    ! -- BCGS WORKING VECTORS
    do n = 1, this%NIABCGS
      this%T(n) = DZERO
      this%V(n) = DZERO
      this%DHAT(n) = DZERO
      this%PHAT(n) = DZERO
      this%QHAT(n) = DZERO
    end do
    !
    ! -- ILUT AND MILUT WORKING VECTORS
    do n = 1, this%NJLU
      this%JLU(n) = IZERO
    end do
    do n = 1, this%NJW
      this%JW(n) = IZERO
    end do
    do n = 1, this%NWLU
      this%WLU(n) = DZERO
    end do
    !
    ! -- REORDERING VECTORS
    do n = 1, i0 + 1
      this%IARO(n) = IZERO
    end do
    do n = 1, iolen
      this%JARO(n) = IZERO
      this%ARO(n) = DZERO
    end do
    !
    ! -- REVERSE CUTHILL MCKEE AND MINIMUM DEGREE ORDERING
    if (this%IORD /= ORD_NONE) then
      call ims_base_calc_order(this%IORD, this%NEQ, this%NJA, this%IA, &
                               this%JA, this%LORDER, this%IORDER)
    end if
    !
    ! -- ALLOCATE MEMORY FOR STORING ITERATION CONVERGENCE DATA
  end subroutine imslinear_ar

  !> @ brief Write summary of settings
    !!
    !!  Write summary of linear accelerator settings.
    !!
  !<
  subroutine imslinear_summary(this, mxiter)
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< ImsLinearDataType instance
    integer(I4B), intent(in) :: mxiter !< maximum number of outer iterations
    ! -- local variables
    character(len=31) :: caccel
    character(len=10) :: clin
    character(len=20) :: cpc
    character(len=20) :: cscale
    character(len=25) :: corder
    character(len=16) :: ccnvg
    character(len=8) :: csmoother
    character(len=15) :: clevel
    character(len=15) :: cdroptol
    integer(I4B) :: i
    integer(I4B) :: j
    ! -- formats
02010 format(1x, /, 7x, 'SOLUTION BY THE', 1x, A31, 1x, 'METHOD', &
           /, 1x, 66('-'), /, &
           ' MAXIMUM OF ', I0, ' CALLS OF SOLUTION ROUTINE', /, &
           ' MAXIMUM OF ', I0, &
           ' INTERNAL ITERATIONS PER CALL TO SOLUTION ROUTINE', /, &
           ' LINEAR ACCELERATION METHOD            =', 1x, A, /, &
           ' MATRIX SCALING APPROACH               =', 1x, A, /, &
           ' MATRIX REORDERING APPROACH            =', 1x, A, /, &
           ' NUMBER OF ORTHOGONALIZATIONS          =', 1x, I0, /, &
           ' HEAD CHANGE CRITERION FOR CLOSURE     =', E15.5, /, &
           ' RESIDUAL CHANGE CRITERION FOR CLOSURE =', E15.5, /, &
           ' RESIDUAL CONVERGENCE NORM             =', 1x, A)
02015 format(' RELAXATION FACTOR (ILU modification)  =', E15.5, /, &
           ' NUMBER OF LEVELS (ILU fill level)     =', A15, /, &
           ' DROP TOLERANCE                        =', A15, //)
02016 format(' RELAXATION FACTOR (cycle correction)  =', E15.5, /, &
           ' NUMBER OF LEVELS (AMG hierarchy)      =', 1x, I0, /, &
           ' AMG SMOOTHER TYPE                     =', 1x, A, /, &
           ' AMG STRENGTH THRESHOLD                =', E15.5, //)
02018 format(/, 7x, A, 1x, 'PRECONDITIONER', &
            /, 1x, 66('-'))
2030 format(1x, A20, 1x, 6(I6, 1x))
2040 format(1x, 20('-'), 1x, 6(6('-'), 1x))
2050 format(1x, 62('-'),/) !
! -- -----------------------------------------------------------
    !
    ! -- resolve integer codes to descriptive strings
    select case (this%ILINMETH)
    case (CG_METHOD)
      caccel = '       CONJUGATE-GRADIENT      '
      clin = 'CG'
    case (BCGS_METHOD)
      caccel = 'BICONJUGATE-GRADIENT STABILIZED'
      clin = 'BCGS'
    case default
      caccel = '             UNKNOWN           '
      clin = 'UNKNOWN'
    end select
    select case (this%IPC)
    case (IPC_ILU0)
      cpc = 'INCOMPLETE LU'
    case (IPC_MILU0)
      cpc = 'MOD. INCOMPLETE LU'
    case (IPC_ILUT)
      cpc = 'INCOMPLETE LUT'
    case (IPC_MILUT)
      cpc = 'MOD. INCOMPLETE LUT'
    case (IPC_AMG)
      cpc = 'AMG'
    case default
      cpc = 'UNKNOWN'
    end select
    select case (this%ISCL)
    case (SCL_DIAGONAL)
      cscale = 'SYMMETRIC SCALING'
    case (SCL_L2NORM)
      cscale = 'L2 NORM SCALING'
    case default
      cscale = 'NO SCALING'
    end select
    select case (this%IORD)
    case (ORD_RCM)
      corder = 'RCM ORDERING'
    case (ORD_MINIMUM_DEGREE)
      corder = 'MINIMUM DEGREE ORDERING'
    case default
      corder = 'ORIGINAL ORDERING'
    end select
    select case (this%ICNVGOPT)
    case (CNVG_INF_NORM_STRICT)
      ccnvg = 'INFINITY NORM S'
    case (CNVG_L2_NORM)
      ccnvg = 'L2 NORM'
    case (CNVG_REL_L2_NORM)
      ccnvg = 'RELATIVE L2 NORM'
    case (CNVG_L2_NORM_REL)
      ccnvg = 'L2 NORM W. REL.'
    case default
      ccnvg = 'INFINITY NORM'
    end select
    select case (this%ISMOOTHER)
    case (SMOOTHER_ILU0)
      csmoother = 'ILU0'
    case default
      csmoother = 'ILU0_ALL'
    end select
    !
    ! -- format clevel and cdroptol (always, so 0 values are shown explicitly)
    write (clevel, '(i15)') this%level
    write (cdroptol, '(e15.5)') this%droptol
    !
    ! -- write solver settings
    write (this%iout, 2010) &
      caccel, MXITER, this%ITER1, &
      trim(clin), trim(cscale), trim(corder), &
      this%NORTH, this%DVCLOSE, this%RCLOSE, &
      trim(ccnvg)
    ! -- write preconditioner-specific settings under a labelled header
    write (this%iout, 2018) trim(cpc)
    if (this%IPC == IPC_AMG) then
      write (this%iout, 2016) this%RELAX, max(1, this%LEVEL), &
        trim(csmoother), this%STHRESH
    else
      write (this%iout, 2015) this%RELAX, trim(adjustl(clevel)), &
        trim(adjustl(cdroptol))
    end if

    if (this%iord /= 0) then
      !
      ! -- WRITE SUMMARY OF REORDERING INFORMATION TO LIST FILE
      if (this%iprims == 2) then
        do i = 1, this%neq, 6
          write (this%iout, 2030) 'ORIGINAL NODE      :', &
            (j, j=i, min(i + 5, this%neq))
          write (this%iout, 2040)
          write (this%iout, 2030) 'REORDERED INDEX    :', &
            (this%lorder(j), j=i, min(i + 5, this%neq))
          write (this%iout, 2030) 'REORDERED NODE     :', &
            (this%iorder(j), j=i, min(i + 5, this%neq))
          write (this%iout, 2050)
        end do
      end if
    end if
  end subroutine imslinear_summary

  !> @ brief Allocate and initialize scalars
    !!
    !!  Allocate and initialize linear accelerator scalars
    !!
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< ImsLinearDataType instance
    !
    ! -- allocate scalars
    call mem_allocate(this%iout, 'IOUT', this%memoryPath)
    call mem_allocate(this%ipc, 'IPC', this%memoryPath)
    call mem_allocate(this%iacpc, 'IACPC', this%memoryPath)
    call mem_allocate(this%niterc, 'NITERC', this%memoryPath)
    call mem_allocate(this%niabcgs, 'NIABCGS', this%memoryPath)
    call mem_allocate(this%niapc, 'NIAPC', this%memoryPath)
    call mem_allocate(this%njapc, 'NJAPC', this%memoryPath)
    call mem_allocate(this%epfact, 'EPFACT', this%memoryPath)
    call mem_allocate(this%l2norm0, 'L2NORM0', this%memoryPath)
    call mem_allocate(this%njlu, 'NJLU', this%memoryPath)
    call mem_allocate(this%njw, 'NJW', this%memoryPath)
    call mem_allocate(this%nwlu, 'NWLU', this%memoryPath)
    !
    ! -- initialize scalars
    this%iout = 0
    this%ipc = 0
    this%iacpc = 0
    this%niterc = 0
    this%niabcgs = 0
    this%niapc = 0
    this%njapc = 0
    this%epfact = DZERO
    this%l2norm0 = 0
    this%njlu = 0
    this%njw = 0
    this%nwlu = 0
  end subroutine allocate_scalars

  !> @ brief Deallocate memory
    !!
    !!  Deallocate linear accelerator memory.
    !!
  !<
  subroutine imslinear_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< linear datatype instance
    !
    !
    ! -- arrays
    call mem_deallocate(this%dscale)
    call mem_deallocate(this%dscale2)
    call precond_destroy(this)
    call mem_deallocate(this%lorder)
    call mem_deallocate(this%iorder)
    call mem_deallocate(this%iaro)
    call mem_deallocate(this%jaro)
    call mem_deallocate(this%aro)
    call mem_deallocate(this%id)
    call mem_deallocate(this%d)
    call mem_deallocate(this%p)
    call mem_deallocate(this%q)
    call mem_deallocate(this%z)
    call mem_deallocate(this%t)
    call mem_deallocate(this%v)
    call mem_deallocate(this%dhat)
    call mem_deallocate(this%phat)
    call mem_deallocate(this%qhat)
    !
    ! -- scalars
    call mem_deallocate(this%iout)
    call mem_deallocate(this%ipc)
    call mem_deallocate(this%iacpc)
    call mem_deallocate(this%niterc)
    call mem_deallocate(this%niabcgs)
    call mem_deallocate(this%niapc)
    call mem_deallocate(this%njapc)
    call mem_deallocate(this%epfact)
    call mem_deallocate(this%l2norm0)
    call mem_deallocate(this%njlu)
    call mem_deallocate(this%njw)
    call mem_deallocate(this%nwlu)
    call mem_deallocate(this%NJA)
    call ims_amg_da(this%amg_prec)
    !
    ! -- nullify pointers
    nullify (this%iprims)
    nullify (this%neq)
    nullify (this%nja)
    nullify (this%ia)
    nullify (this%ja)
    nullify (this%amat)
    nullify (this%rhs)
    nullify (this%x)
  end subroutine imslinear_da

  !> @brief Allocate and initialize the preconditioner work arrays
  !!
  !! Determines the preconditioner array dimensions from the current
  !! NEQ/NJA/LEVEL/IPC and allocates the ILU-family work arrays
  !! (IAPC/JAPC/APC/IW/W/JLU/JW/WLU); for ILU0/MILU0 the IAPC/JAPC sparsity is
  !! generated. Separated from imslinear_ar so the preconditioner can be
  !! (re)allocated at runtime.
  !<
  subroutine precond_allocate(this)
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< ImsLinearDataType instance
    ! -- local variables
    integer(I4B) :: n
    !
    ! -- release any existing preconditioner arrays so this routine is
    !    idempotent and can be called again to reallocate at runtime
    if (associated(this%IAPC)) call precond_destroy(this)
    !
    ! -- determine dimensions for preconditioning arrays
    call ims_calc_pcdims(this%NEQ, this%NJA, this%IA, this%LEVEL, this%IPC, &
                         this%NIAPC, this%NJAPC, this%NJLU, this%NJW, this%NWLU)
    !
    ! -- allocate base preconditioner vectors
    call mem_allocate(this%IAPC, this%NIAPC + 1, 'IAPC', trim(this%memoryPath))
    call mem_allocate(this%JAPC, this%NJAPC, 'JAPC', trim(this%memoryPath))
    call mem_allocate(this%APC, this%NJAPC, 'APC', trim(this%memoryPath))
    !
    ! -- allocate ILU0/MILU0 non-zero row entry vectors
    call mem_allocate(this%IW, this%NIAPC, 'IW', trim(this%memoryPath))
    call mem_allocate(this%W, this%NIAPC, 'W', trim(this%memoryPath))
    !
    ! -- allocate ILUT vectors
    call mem_allocate(this%JLU, this%NJLU, 'JLU', trim(this%memoryPath))
    call mem_allocate(this%JW, this%NJW, 'JW', trim(this%memoryPath))
    call mem_allocate(this%WLU, this%NWLU, 'WLU', trim(this%memoryPath))
    !
    ! -- generate IAPC and JAPC for ILU0 and MILU0
    if (this%IPC == IPC_ILU0 .or. this%IPC == IPC_MILU0) then
      call ims_base_pccrs(this%NEQ, this%NJA, this%IA, this%JA, &
                          this%IAPC, this%JAPC)
    end if
    !
    ! -- initialize preconditioner values
    do n = 1, this%NJAPC
      this%APC(n) = DZERO
    end do
  end subroutine precond_allocate

  !> @brief Deallocate the preconditioner work arrays
  !!
  !! Preconditioner-only teardown, separated from imslinear_da so the
  !! preconditioner can be destroyed and reallocated at runtime without
  !! disturbing the solver work arrays.
  !<
  subroutine precond_destroy(this)
    use MemoryManagerExtModule, only: memorystore_release
    ! -- dummy variables
    class(ImsLinearDataType), intent(inout) :: this !< ImsLinearDataType instance
    !
    ! -- mem_deallocate is a no-op in the memory manager, so releasing the named
    !    store entries (and nullifying the aliases below) is what actually frees
    !    these arrays and lets the preconditioner be reallocated at runtime
    call memorystore_release('IAPC', trim(this%memoryPath))
    call memorystore_release('JAPC', trim(this%memoryPath))
    call memorystore_release('APC', trim(this%memoryPath))
    call memorystore_release('IW', trim(this%memoryPath))
    call memorystore_release('W', trim(this%memoryPath))
    call memorystore_release('JLU', trim(this%memoryPath))
    call memorystore_release('JW', trim(this%memoryPath))
    call memorystore_release('WLU', trim(this%memoryPath))
    !
    nullify (this%IAPC)
    nullify (this%JAPC)
    nullify (this%APC)
    nullify (this%IW)
    nullify (this%W)
    nullify (this%JLU)
    nullify (this%JW)
    nullify (this%WLU)
  end subroutine precond_destroy

  !> @ brief Set default settings
    !!
    !!  Set default linear accelerator settings.
    !!
  !<
  subroutine imslinear_set_input(this, IFDPARAM)
    ! -- dummy variables
    class(ImsLinearDataType), intent(INOUT) :: this !< ImsLinearDataType instance
    integer(I4B), intent(IN) :: IFDPARAM !< complexity option
    ! -- code
    select case (IFDPARAM)
      !
      ! -- Simple option
    case (1)
      this%ITER1 = 50
      this%ILINMETH = 1
      this%IPC = IPC_ILU0
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM3
      this%RCLOSE = DEM1
      this%RELAX = DZERO
      this%LEVEL = 0
      this%DROPTOL = DZERO
      this%NORTH = 0
      !
      ! -- Moderate
    case (2)
      this%ITER1 = 100
      this%ILINMETH = 2
      this%IPC = IPC_MILU0
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM2
      this%RCLOSE = DEM1
      this%RELAX = 0.97d0
      this%LEVEL = 0
      this%DROPTOL = DZERO
      this%NORTH = 0
      !
      ! -- Complex
    case (3)
      this%ITER1 = 500
      this%ILINMETH = 2
      this%IPC = IPC_ILUT
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM1
      this%RCLOSE = DEM1
      this%RELAX = DZERO
      this%LEVEL = 5
      this%DROPTOL = DEM4
      this%NORTH = 2
    end select
  end subroutine imslinear_set_input

  !> @ brief Base linear accelerator subroutine
    !!
    !!  Base linear accelerator subroutine that scales and reorders
    !!  the system of equations, if necessary, updates the preconditioner,
    !!  and calls the appropriate linear accelerator.
    !!
  !<
  subroutine imslinear_ap(this, ICNVG, KSTP, KITER, IN_ITER, &
                          NCONV, CONVNMOD, CONVMODSTART, &
                          CACCEL, summary)
    ! -- modules
    use SimModule
    ! -- dummy variables
    class(ImsLinearDataType), intent(INOUT) :: this !< ImsLinearDataType instance
    integer(I4B), intent(INOUT) :: ICNVG !< convergence flag (1) non-convergence (0)
    integer(I4B), intent(IN) :: KSTP !< time step number
    integer(I4B), intent(IN) :: KITER !< outer iteration number
    integer(I4B), intent(INOUT) :: IN_ITER !< inner iteration number
    ! -- convergence information dummy variables
    integer(I4B), intent(IN) :: NCONV !<
    integer(I4B), intent(IN) :: CONVNMOD !<
    integer(I4B), dimension(CONVNMOD + 1), intent(INOUT) :: CONVMODSTART !<
    character(len=31), dimension(NCONV), intent(INOUT) :: CACCEL !<
    type(ConvergenceSummaryType), pointer, intent(in) :: summary !< Convergence summary report
    ! -- local variables
    integer(I4B) :: n
    integer(I4B) :: innerit
    integer(I4B) :: irc
    integer(I4B) :: itmax
    integer(I4B) :: nlev_prev
    real(DP) :: dnrm2
    !
    ! -- set epfact based on timestep
    this%EPFACT = ims_base_epfact(this%ICNVGOPT, KSTP)
    !
    ! -- SCALE PROBLEM
    if (this%ISCL /= SCL_NONE) then
      call ims_base_scale(0, this%ISCL, &
                          this%NEQ, this%NJA, this%IA, this%JA, &
                          this%AMAT, this%X, this%RHS, &
                          this%DSCALE, this%DSCALE2)
    end if
    !
    ! -- PERMUTE ROWS, COLUMNS, AND RHS
    if (this%IORD /= ORD_NONE) then
      call dperm(this%NEQ, this%AMAT, this%JA, this%IA, &
                 this%ARO, this%JARO, this%IARO, &
                 this%LORDER, this%ID, 1)
      call dvperm(this%NEQ, this%X, this%LORDER)
      call dvperm(this%NEQ, this%RHS, this%LORDER)
      this%IA0 => this%IARO
      this%JA0 => this%JARO
      this%A0 => this%ARO
    else
      this%IA0 => this%IA
      this%JA0 => this%JA
      this%A0 => this%AMAT
    end if
    !
    ! -- UPDATE PRECONDITIONER
    if (this%IPC == IPC_AMG) then
      call g_prof%start("AMG setup", this%amg_prec%itmr_setup)
      nlev_prev = this%amg_prec%nlevels
      call ims_amg_setup(this%amg_prec, this%NEQ, this%NJA, this%IA0, this%JA0, &
                         this%A0, max(1, this%LEVEL), this%NSMOOTH, &
                         this%RELAX, this%STHRESH, this%ISMOOTHER)
      call g_prof%stop(this%amg_prec%itmr_setup)
      ! report the hierarchy the first time it is built; later calls reuse the
      ! structure and only update the matrix values
      if (nlev_prev == 0) then
        call ims_amg_summary(this%amg_prec, this%iout, this%STHRESH)
      end if
    else
      call ims_base_pcu(this%iout, this%NJA, this%NEQ, this%NIAPC, this%NJAPC, &
                        this%IPC, this%RELAX, this%A0, this%IA0, this%JA0, &
                        this%APC, this%IAPC, this%JAPC, this%IW, this%W, &
                        this%LEVEL, this%DROPTOL, this%NJLU, this%NJW, &
                        this%NWLU, this%JLU, this%JW, this%WLU)
    end if
    !
    ! -- INITIALIZE SOLUTION VARIABLE AND ARRAYS
    if (KITER == 1) then
      this%NITERC = 0
      summary%iter_cnt = 0
    end if
    irc = 1
    ICNVG = 0
    do n = 1, this%NEQ
      this%D(n) = DZERO
      this%P(n) = DZERO
      this%Q(n) = DZERO
      this%Z(n) = DZERO
    end do
    !
    ! -- CALCULATE INITIAL RESIDUAL
    call ims_base_residual(this%NEQ, this%NJA, this%X, this%RHS, this%D, &
                           this%A0, this%IA0, this%JA0)
    this%L2NORM0 = dnrm2(this%NEQ, this%D, 1)
    !
    ! -- CHECK FOR EXACT SOLUTION
    itmax = this%ITER1
    if (this%L2NORM0 == DZERO) then
      itmax = 0
      ICNVG = 1
    end if
    !
    ! -- SOLUTION BY THE CONJUGATE GRADIENT METHOD
    if (this%ILINMETH == CG_METHOD) then
      call ims_base_cg(ICNVG, itmax, innerit, &
                       this%NEQ, this%NJA, this%NIAPC, this%NJAPC, &
                       this%IPC, this%ICNVGOPT, this%NORTH, &
                       this%DVCLOSE, this%RCLOSE, this%L2NORM0, &
                       this%EPFACT, this%IA0, this%JA0, this%A0, &
                       this%IAPC, this%JAPC, this%APC, &
                       this%X, this%RHS, this%D, this%P, this%Q, this%Z, &
                       this%NJLU, this%IW, this%JLU, &
                       NCONV, CONVNMOD, CONVMODSTART, &
                       CACCEL, summary, amg=this%amg_prec)
      !
      ! -- SOLUTION BY THE BICONJUGATE GRADIENT STABILIZED METHOD
    else if (this%ILINMETH == BCGS_METHOD) then
      call ims_base_bcgs(ICNVG, itmax, innerit, &
                         this%NEQ, this%NJA, this%NIAPC, this%NJAPC, &
                         this%IPC, this%ICNVGOPT, this%NORTH, &
                         this%ISCL, this%DSCALE, &
                         this%DVCLOSE, this%RCLOSE, this%L2NORM0, &
                         this%EPFACT, this%IA0, this%JA0, this%A0, &
                         this%IAPC, this%JAPC, this%APC, &
                         this%X, this%RHS, this%D, this%P, this%Q, &
                         this%T, this%V, this%DHAT, this%PHAT, this%QHAT, &
                         this%NJLU, this%IW, this%JLU, &
                         NCONV, CONVNMOD, CONVMODSTART, &
                         CACCEL, summary, amg=this%amg_prec)
    end if
    !
    ! -- BACK PERMUTE AMAT, SOLUTION, AND RHS
    if (this%IORD /= ORD_NONE) then
      call dperm(this%NEQ, this%A0, this%JA0, this%IA0, &
                 this%AMAT, this%JA, this%IA, &
                 this%IORDER, this%ID, 1)
      call dvperm(this%NEQ, this%X, this%IORDER)
      call dvperm(this%NEQ, this%RHS, this%IORDER)
    end if
    !
    ! -- UNSCALE PROBLEM
    if (this%ISCL /= SCL_NONE) then
      call ims_base_scale(1, this%ISCL, &
                          this%NEQ, this%NJA, this%IA, this%JA, &
                          this%AMAT, this%X, this%RHS, &
                          this%DSCALE, this%DSCALE2)
    end if
    !
    ! -- SET IMS INNER ITERATION NUMBER (IN_ITER) TO NUMBER OF
    !       IMSLINEAR INNER ITERATIONS (innerit)
    IN_ITER = innerit
  end subroutine imslinear_ap

end module IMSLinearModule
