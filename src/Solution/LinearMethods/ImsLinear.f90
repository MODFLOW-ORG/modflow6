MODULE IMSLinearModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, LENSOLUTIONNAME, LENMEMPATH, &
                             IZERO, DZERO, DPREC, DSAME, &
                             DEM8, DEM6, DEM5, DEM4, DEM3, DEM2, DEM1, &
                             DHALF, DONE, DTWO, &
                             VDEBUG
  use IMSLinearBaseModule, only: ims_base_cg, ims_base_bcgs, &
                                 ims_base_pccrs, ims_base_calc_order, &
                                 ims_base_scale, ims_base_pcu, &
                                 ims_base_residual, ims_base_epfact, &
                                 ims_calc_pcdims
  use BlockParserModule, only: BlockParserType
  use MatrixBaseModule
  use ConvergenceSummaryModule
  use ImsLinearSettingsModule

  IMPLICIT NONE
  private

  TYPE, PUBLIC :: ImsLinearDataType
    character(len=LENMEMPATH) :: memoryPath !< the path for storing variables in the memory manager
    integer(I4B), POINTER :: iout => NULL() !< simulation listing file unit
    integer(I4B), POINTER :: IPRIMS => NULL() !< print flag
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
    !
    integer(I4B), POINTER :: IPC => NULL() !< preconditioner flag
    integer(I4B), POINTER :: IACPC => NULL() !< preconditioner CRS row pointers
    integer(I4B), POINTER :: NITERC => NULL() !<
    integer(I4B), POINTER :: NIABCGS => NULL() !< size of working vectors for BCGS linear accelerator
    integer(I4B), POINTER :: NIAPC => NULL() !< preconditioner number of rows
    integer(I4B), POINTER :: NJAPC => NULL() !< preconditioner number of non-zero entries
    real(DP), POINTER :: EPFACT => NULL() !< factor for decreasing convergence criteria in seubsequent Picard iterations
    real(DP), POINTER :: L2NORM0 => NULL() !< initial L2 norm
    ! -- ilut variables
    integer(I4B), POINTER :: NJLU => NULL() !< length of jlu work vector
    integer(I4B), POINTER :: NJW => NULL() !< length of jw work vector
    integer(I4B), POINTER :: NWLU => NULL() !< length of wlu work vector
    ! -- pointers to solution variables
    integer(I4B), POINTER :: NEQ => NULL() !< number of equations (rows in matrix)
    integer(I4B), POINTER :: NJA => NULL() !< number of non-zero values in amat
    integer(I4B), dimension(:), pointer, contiguous :: IA => NULL() !< position of start of each row
    integer(I4B), dimension(:), pointer, contiguous :: JA => NULL() !< column pointer
    real(DP), dimension(:), pointer, contiguous :: AMAT => NULL() !< coefficient matrix
    real(DP), dimension(:), pointer, contiguous :: RHS => NULL() !< right-hand side of equation
    real(DP), dimension(:), pointer, contiguous :: X => NULL() !< dependent variable
    ! VECTORS
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: DSCALE => NULL() !< scaling factor
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: DSCALE2 => NULL() !< unscaling factor
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IAPC => NULL() !< position of start of each row in preconditioner matrix
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JAPC => NULL() !< preconditioner matrix column pointer
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: APC => NULL() !< preconditioner coefficient matrix
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: LORDER => NULL() !< reordering mapping
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IORDER => NULL() !< mapping to restore reordered matrix
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IARO => NULL() !< position of start of each row in reordered matrix
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JARO => NULL() !< reordered matrix column pointer
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: ARO => NULL() !< reordered coefficient matrix
    ! WORKING ARRAYS
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IW => NULL() !< integer working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: W => NULL() !< real working array
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: ID => NULL() !< integer working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: D => NULL() !< real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: P => NULL() !< real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: Q => NULL() !< real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: Z => NULL() !< real working array
    ! BICGSTAB WORKING ARRAYS
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: T => NULL() !< BICGSTAB real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: V => NULL() !< BICGSTAB real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: DHAT => NULL() !< BICGSTAB real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: PHAT => NULL() !< BICGSTAB real working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: QHAT => NULL() !< rBICGSTAB eal working array
    ! POINTERS FOR USE WITH BOTH ORIGINAL AND RCM ORDERINGS
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: IA0 => NULL() !< pointer to current CRS row pointers
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JA0 => NULL() !< pointer to current CRS column pointers
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: A0 => NULL() !< pointer to current coefficient matrix
    ! ILUT WORKING ARRAYS
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JLU => NULL() !< ilut integer working array
    integer(I4B), POINTER, DIMENSION(:), CONTIGUOUS :: JW => NULL() !< ilut integer working array
    real(DP), POINTER, DIMENSION(:), CONTIGUOUS :: WLU => NULL() !< ilut real working array

    ! PROCEDURES (METHODS)
  CONTAINS
    PROCEDURE :: IMSLINEAR_ALLOCATE => imslinear_ar
    procedure :: imslinear_summary
    PROCEDURE :: IMSLINEAR_APPLY => imslinear_ap
    procedure :: IMSLINEAR_DA => imslinear_da
    procedure, private :: allocate_scalars
    ! -- PRIVATE PROCEDURES
    PROCEDURE, PRIVATE :: SET_IMSLINEAR_INPUT => imslinear_set_input
  END TYPE ImsLinearDataType

CONTAINS

  !> @ brief Allocate storage and read data
    !!
    !!  Allocate storage for linear accelerators and read data
    !!
  !<
  SUBROUTINE imslinear_ar(this, NAME, IOUT, IPRIMS, MXITER, &
                          NEQ, matrix, RHS, X, linear_settings)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use MemoryHelperModule, only: create_mem_path
    use SimModule, only: store_error, count_errors, &
                         deprecation_warning
    ! -- dummy variables
    CLASS(ImsLinearDataType), INTENT(INOUT) :: this !< ImsLinearDataType instance
    CHARACTER(LEN=LENSOLUTIONNAME), INTENT(IN) :: NAME !< solution name
    integer(I4B), INTENT(IN) :: IOUT !< simulation listing file unit
    integer(I4B), TARGET, INTENT(IN) :: IPRIMS !< print option
    integer(I4B), INTENT(IN) :: MXITER !< maximum outer iterations
    integer(I4B), TARGET, INTENT(IN) :: NEQ !< number of equations
    class(MatrixBaseType), pointer :: matrix
    real(DP), DIMENSION(NEQ), TARGET, INTENT(INOUT) :: RHS !< right-hand side
    real(DP), DIMENSION(NEQ), TARGET, INTENT(INOUT) :: X !< dependent variables
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
    this%IPC = 0
    !
    this%IACPC = 0
    !
    ! -- PRINT A MESSAGE IDENTIFYING IMSLINEAR SOLVER PACKAGE
    write (iout, 2000)
02000 FORMAT(1X, /1X, 'IMSLINEAR -- UNSTRUCTURED LINEAR SOLUTION', &
           ' PACKAGE, VERSION 8, 04/28/2017')
    !
    ! -- DETERMINE PRECONDITIONER
    IF (this%LEVEL > 0 .OR. this%DROPTOL > DZERO) THEN
      this%IPC = 3
    ELSE
      this%IPC = 1
    END IF
    IF (this%RELAX > DZERO) THEN
      this%IPC = this%IPC + 1
    END IF
    !
    ! -- ERROR CHECKING FOR OPTIONS
    IF (this%ISCL < 0) this%ISCL = 0
    IF (this%ISCL > 2) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR ISCL MUST BE <= 2'
      call store_error(errmsg)
    END IF
    IF (this%IORD < 0) this%IORD = 0
    IF (this%IORD > 2) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR IORD MUST BE <= 2'
      call store_error(errmsg)
    END IF
    IF (this%NORTH < 0) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR NORTH MUST >= 0'
      call store_error(errmsg)
    END IF
    IF (this%RCLOSE == DZERO) THEN
      IF (this%ICNVGOPT /= 3) THEN
        WRITE (errmsg, '(A)') 'IMSLINEAR7AR RCLOSE MUST > 0.0'
        call store_error(errmsg)
      END IF
    END IF
    IF (this%RELAX < DZERO) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR RELAX MUST BE >= 0.0'
      call store_error(errmsg)
    END IF
    IF (this%RELAX > DONE) THEN
      WRITE (errmsg, '(A)') 'IMSLINEAR7AR RELAX MUST BE <= 1.0'
      call store_error(errmsg)
    END IF
    !
    ! -- INITIALIZE IMSLINEAR VARIABLES
    this%NITERC = 0
    !
    ! -- ALLOCATE AND INITIALIZE MEMORY FOR IMSLINEAR
    iscllen = 1
    IF (this%ISCL .NE. 0) iscllen = NEQ
    CALL mem_allocate(this%DSCALE, iscllen, 'DSCALE', TRIM(this%memoryPath))
    CALL mem_allocate(this%DSCALE2, iscllen, 'DSCALE2', TRIM(this%memoryPath))
    !
    ! -- determine dimensions for preconditing arrays
    call ims_calc_pcdims(this%NEQ, this%NJA, this%IA, this%LEVEL, this%IPC, &
                         this%NIAPC, this%NJAPC, this%NJLU, this%NJW, this%NWLU)
    !
    ! -- ALLOCATE BASE PRECONDITIONER VECTORS
    CALL mem_allocate(this%IAPC, this%NIAPC + 1, 'IAPC', TRIM(this%memoryPath))
    CALL mem_allocate(this%JAPC, this%NJAPC, 'JAPC', TRIM(this%memoryPath))
    CALL mem_allocate(this%APC, this%NJAPC, 'APC', TRIM(this%memoryPath))
    !
    ! -- ALLOCATE MEMORY FOR ILU0 AND MILU0 NON-ZERO ROW ENTRY VECTOR
    CALL mem_allocate(this%IW, this%NIAPC, 'IW', TRIM(this%memoryPath))
    CALL mem_allocate(this%W, this%NIAPC, 'W', TRIM(this%memoryPath))
    !
    ! -- ALLOCATE MEMORY FOR ILUT VECTORS
    CALL mem_allocate(this%JLU, this%NJLU, 'JLU', TRIM(this%memoryPath))
    CALL mem_allocate(this%JW, this%NJW, 'JW', TRIM(this%memoryPath))
    CALL mem_allocate(this%WLU, this%NWLU, 'WLU', TRIM(this%memoryPath))
    !
    ! -- GENERATE IAPC AND JAPC FOR ILU0 AND MILU0
    IF (this%IPC == 1 .OR. this%IPC == 2) THEN
      CALL ims_base_pccrs(this%NEQ, this%NJA, this%IA, this%JA, &
                          this%IAPC, this%JAPC)
    END IF
    !
    ! -- ALLOCATE SPACE FOR PERMUTATION VECTOR
    i0 = 1
    iolen = 1
    IF (this%IORD .NE. 0) THEN
      i0 = this%NEQ
      iolen = this%NJA
    END IF
    CALL mem_allocate(this%LORDER, i0, 'LORDER', TRIM(this%memoryPath))
    CALL mem_allocate(this%IORDER, i0, 'IORDER', TRIM(this%memoryPath))
    CALL mem_allocate(this%IARO, i0 + 1, 'IARO', TRIM(this%memoryPath))
    CALL mem_allocate(this%JARO, iolen, 'JARO', TRIM(this%memoryPath))
    CALL mem_allocate(this%ARO, iolen, 'ARO', TRIM(this%memoryPath))
    !
    ! -- ALLOCATE WORKING VECTORS FOR IMSLINEAR SOLVER
    CALL mem_allocate(this%ID, this%NEQ, 'ID', TRIM(this%memoryPath))
    CALL mem_allocate(this%D, this%NEQ, 'D', TRIM(this%memoryPath))
    CALL mem_allocate(this%P, this%NEQ, 'P', TRIM(this%memoryPath))
    CALL mem_allocate(this%Q, this%NEQ, 'Q', TRIM(this%memoryPath))
    CALL mem_allocate(this%Z, this%NEQ, 'Z', TRIM(this%memoryPath))
    !
    ! -- ALLOCATE MEMORY FOR BCGS WORKING ARRAYS
    this%NIABCGS = 1
    IF (this%ILINMETH == 2) THEN
      this%NIABCGS = this%NEQ
    END IF
    CALL mem_allocate(this%T, this%NIABCGS, 'T', TRIM(this%memoryPath))
    CALL mem_allocate(this%V, this%NIABCGS, 'V', TRIM(this%memoryPath))
    CALL mem_allocate(this%DHAT, this%NIABCGS, 'DHAT', TRIM(this%memoryPath))
    CALL mem_allocate(this%PHAT, this%NIABCGS, 'PHAT', TRIM(this%memoryPath))
    CALL mem_allocate(this%QHAT, this%NIABCGS, 'QHAT', TRIM(this%memoryPath))
    !
    ! -- INITIALIZE IMSLINEAR VECTORS
    DO n = 1, iscllen
      this%DSCALE(n) = DONE
      this%DSCALE2(n) = DONE
    END DO
    DO n = 1, this%NJAPC
      this%APC(n) = DZERO
    END DO
    !
    ! -- WORKING VECTORS
    DO n = 1, this%NEQ
      this%ID(n) = IZERO
      this%D(n) = DZERO
      this%P(n) = DZERO
      this%Q(n) = DZERO
      this%Z(n) = DZERO
    END DO
    DO n = 1, this%NIAPC
      this%IW(n) = IZERO
      this%W(n) = DZERO
    END DO
    !
    ! -- BCGS WORKING VECTORS
    DO n = 1, this%NIABCGS
      this%T(n) = DZERO
      this%V(n) = DZERO
      this%DHAT(n) = DZERO
      this%PHAT(n) = DZERO
      this%QHAT(n) = DZERO
    END DO
    !
    ! -- ILUT AND MILUT WORKING VECTORS
    DO n = 1, this%NJLU
      this%JLU(n) = IZERO
    END DO
    DO n = 1, this%NJW
      this%JW(n) = IZERO
    END DO
    DO n = 1, this%NWLU
      this%WLU(n) = DZERO
    END DO
    !
    ! -- REORDERING VECTORS
    DO n = 1, i0 + 1
      this%IARO(n) = IZERO
    END DO
    DO n = 1, iolen
      this%JARO(n) = IZERO
      this%ARO(n) = DZERO
    END DO
    !
    ! -- REVERSE CUTHILL MCKEE AND MINIMUM DEGREE ORDERING
    IF (this%IORD .NE. 0) THEN
      CALL ims_base_calc_order(this%IORD, this%NEQ, this%NJA, this%IA, &
                               this%JA, this%LORDER, this%IORDER)
    END IF
    !
    ! -- ALLOCATE MEMORY FOR STORING ITERATION CONVERGENCE DATA
  end SUBROUTINE imslinear_ar

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
    CHARACTER(LEN=10) :: clin(0:2)
    CHARACTER(LEN=31) :: clintit(0:2)
    CHARACTER(LEN=20) :: cipc(0:4)
    CHARACTER(LEN=20) :: cscale(0:2)
    CHARACTER(LEN=25) :: corder(0:2)
    CHARACTER(LEN=16), DIMENSION(0:4) :: ccnvgopt
    CHARACTER(LEN=15) :: clevel
    CHARACTER(LEN=15) :: cdroptol
    integer(I4B) :: i
    integer(I4B) :: j
    ! -- data
    DATA clin/'UNKNOWN   ', &
             &'CG        ', &
             &'BCGS      '/
    DATA clintit/'             UNKNOWN           ', &
                &'       CONJUGATE-GRADIENT      ', &
                &'BICONJUGATE-GRADIENT STABILIZED'/
    DATA cipc/'UNKNOWN             ', &
             &'INCOMPLETE LU       ', &
             &'MOD. INCOMPLETE LU  ', &
             &'INCOMPLETE LUT      ', &
             &'MOD. INCOMPLETE LUT '/
    DATA cscale/'NO SCALING          ', &
               &'SYMMETRIC SCALING   ', &
               &'L2 NORM SCALING     '/
    DATA corder/'ORIGINAL ORDERING        ', &
               &'RCM ORDERING             ', &
               &'MINIMUM DEGREE ORDERING  '/
    DATA ccnvgopt/'INFINITY NORM   ', &
                 &'INFINITY NORM S ', &
                 &'L2 NORM         ', &
                 &'RELATIVE L2NORM ', &
                 &'L2 NORM W. REL. '/
    ! -- formats
02010 FORMAT(1X, /, 7X, 'SOLUTION BY THE', 1X, A31, 1X, 'METHOD', &
           /, 1X, 66('-'), /, &
           ' MAXIMUM OF ', I0, ' CALLS OF SOLUTION ROUTINE', /, &
           ' MAXIMUM OF ', I0, &
           ' INTERNAL ITERATIONS PER CALL TO SOLUTION ROUTINE', /, &
           ' LINEAR ACCELERATION METHOD            =', 1X, A, /, &
           ' MATRIX PRECONDITIONING TYPE           =', 1X, A, /, &
           ' MATRIX SCALING APPROACH               =', 1X, A, /, &
           ' MATRIX REORDERING APPROACH            =', 1X, A, /, &
           ' NUMBER OF ORTHOGONALIZATIONS          =', 1X, I0, /, &
           ' HEAD CHANGE CRITERION FOR CLOSURE     =', E15.5, /, &
           ' RESIDUAL CHANGE CRITERION FOR CLOSURE =', E15.5, /, &
           ' RESIDUAL CONVERGENCE OPTION           =', 1X, I0, /, &
           ' RESIDUAL CONVERGENCE NORM             =', 1X, A, /, &
           ' RELAXATION FACTOR                     =', E15.5)
02015 FORMAT(' NUMBER OF LEVELS                      =', A15, /, &
           ' DROP TOLERANCE                        =', A15, //)
2030 FORMAT(1X, A20, 1X, 6(I6, 1X))
2040 FORMAT(1X, 20('-'), 1X, 6(6('-'), 1X))
2050 FORMAT(1X, 62('-'),/) !
! -- -----------------------------------------------------------
    !
    ! -- initialize clevel and cdroptol
    clevel = ''
    cdroptol = ''
    !
    ! -- write common variables to all linear accelerators
    write (this%iout, 2010) &
      clintit(this%ILINMETH), MXITER, this%ITER1, &
      clin(this%ILINMETH), cipc(this%IPC), &
      cscale(this%ISCL), corder(this%IORD), &
      this%NORTH, this%DVCLOSE, this%RCLOSE, &
      this%ICNVGOPT, ccnvgopt(this%ICNVGOPT), &
      this%RELAX
    if (this%level > 0) then
      write (clevel, '(i15)') this%level
    end if
    if (this%droptol > DZERO) then
      write (cdroptol, '(e15.5)') this%droptol
    end if
    IF (this%level > 0 .or. this%droptol > DZERO) THEN
      write (this%iout, 2015) trim(adjustl(clevel)), &
        trim(adjustl(cdroptol))
    ELSE
      write (this%iout, '(//)')
    END IF

    if (this%iord /= 0) then
      !
      ! -- WRITE SUMMARY OF REORDERING INFORMATION TO LIST FILE
      if (this%iprims == 2) then
        DO i = 1, this%neq, 6
          write (this%iout, 2030) 'ORIGINAL NODE      :', &
            (j, j=i, MIN(i + 5, this%neq))
          write (this%iout, 2040)
          write (this%iout, 2030) 'REORDERED INDEX    :', &
            (this%lorder(j), j=i, MIN(i + 5, this%neq))
          write (this%iout, 2030) 'REORDERED NODE     :', &
            (this%iorder(j), j=i, MIN(i + 5, this%neq))
          write (this%iout, 2050)
        END DO
      END IF
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
    ! -- arrays
    call mem_deallocate(this%dscale)
    call mem_deallocate(this%dscale2)
    call mem_deallocate(this%iapc)
    call mem_deallocate(this%japc)
    call mem_deallocate(this%apc)
    call mem_deallocate(this%iw)
    call mem_deallocate(this%w)
    call mem_deallocate(this%jlu)
    call mem_deallocate(this%jw)
    call mem_deallocate(this%wlu)
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

  !> @ brief Set default settings
    !!
    !!  Set default linear accelerator settings.
    !!
  !<
  SUBROUTINE imslinear_set_input(this, IFDPARAM)
    ! -- dummy variables
    CLASS(ImsLinearDataType), INTENT(INOUT) :: this !< ImsLinearDataType instance
    integer(I4B), INTENT(IN) :: IFDPARAM !< complexity option
    ! -- code
    SELECT CASE (IFDPARAM)
      !
      ! -- Simple option
    CASE (1)
      this%ITER1 = 50
      this%ILINMETH = 1
      this%IPC = 1
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
    CASE (2)
      this%ITER1 = 100
      this%ILINMETH = 2
      this%IPC = 2
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM2
      this%RCLOSE = DEM1
      this%RELAX = 0.97D0
      this%LEVEL = 0
      this%DROPTOL = DZERO
      this%NORTH = 0
      !
      ! -- Complex
    CASE (3)
      this%ITER1 = 500
      this%ILINMETH = 2
      this%IPC = 3
      this%ISCL = 0
      this%IORD = 0
      this%DVCLOSE = DEM1
      this%RCLOSE = DEM1
      this%RELAX = DZERO
      this%LEVEL = 5
      this%DROPTOL = DEM4
      this%NORTH = 2
    END SELECT
  end SUBROUTINE imslinear_set_input

  !> @ brief Base linear accelerator subroutine
    !!
    !!  Base linear accelerator subroutine that scales and reorders
    !!  the system of equations, if necessary, updates the preconditioner,
    !!  and calls the appropriate linear accelerator.
    !!
  !<
  SUBROUTINE imslinear_ap(this, ICNVG, KSTP, KITER, IN_ITER, &
                          NCONV, CONVNMOD, CONVMODSTART, &
                          CACCEL, summary)
    ! -- modules
    USE SimModule
    ! -- dummy variables
    CLASS(ImsLinearDataType), INTENT(INOUT) :: this !< ImsLinearDataType instance
    integer(I4B), INTENT(INOUT) :: ICNVG !< convergence flag (1) non-convergence (0)
    integer(I4B), INTENT(IN) :: KSTP !< time step number
    integer(I4B), INTENT(IN) :: KITER !< outer iteration number
    integer(I4B), INTENT(INOUT) :: IN_ITER !< inner iteration number
    ! -- convergence information dummy variables
    integer(I4B), INTENT(IN) :: NCONV !<
    integer(I4B), INTENT(IN) :: CONVNMOD !<
    integer(I4B), DIMENSION(CONVNMOD + 1), INTENT(INOUT) :: CONVMODSTART !<
    character(len=31), DIMENSION(NCONV), INTENT(INOUT) :: CACCEL !<
    type(ConvergenceSummaryType), pointer, intent(in) :: summary !< Convergence summary report
    ! -- local variables
    integer(I4B) :: n
    integer(I4B) :: innerit
    integer(I4B) :: irc
    integer(I4B) :: itmax
    real(DP) :: dnrm2
    !
    ! -- set epfact based on timestep
    this%EPFACT = ims_base_epfact(this%ICNVGOPT, KSTP)
    !
    ! -- SCALE PROBLEM
    IF (this%ISCL .NE. 0) THEN
      CALL ims_base_scale(0, this%ISCL, &
                          this%NEQ, this%NJA, this%IA, this%JA, &
                          this%AMAT, this%X, this%RHS, &
                          this%DSCALE, this%DSCALE2)
    END IF
    !
    ! -- PERMUTE ROWS, COLUMNS, AND RHS
    IF (this%IORD /= 0) THEN
      CALL dperm(this%NEQ, this%AMAT, this%JA, this%IA, &
                 this%ARO, this%JARO, this%IARO, &
                 this%LORDER, this%ID, 1)
      CALL dvperm(this%NEQ, this%X, this%LORDER)
      CALL dvperm(this%NEQ, this%RHS, this%LORDER)
      this%IA0 => this%IARO
      this%JA0 => this%JARO
      this%A0 => this%ARO
    ELSE
      this%IA0 => this%IA
      this%JA0 => this%JA
      this%A0 => this%AMAT
    END IF
    !
    ! -- UPDATE PRECONDITIONER
    CALL ims_base_pcu(this%iout, this%NJA, this%NEQ, this%NIAPC, this%NJAPC, &
                      this%IPC, this%RELAX, this%A0, this%IA0, this%JA0, &
                      this%APC, this%IAPC, this%JAPC, this%IW, this%W, &
                      this%LEVEL, this%DROPTOL, this%NJLU, this%NJW, &
                      this%NWLU, this%JLU, this%JW, this%WLU)
    !
    ! -- INITIALIZE SOLUTION VARIABLE AND ARRAYS
    IF (KITER == 1) then
      this%NITERC = 0
      summary%iter_cnt = 0
    end if
    irc = 1
    ICNVG = 0
    DO n = 1, this%NEQ
      this%D(n) = DZERO
      this%P(n) = DZERO
      this%Q(n) = DZERO
      this%Z(n) = DZERO
    END DO
    !
    ! -- CALCULATE INITIAL RESIDUAL
    call ims_base_residual(this%NEQ, this%NJA, this%X, this%RHS, this%D, &
                           this%A0, this%IA0, this%JA0)
    this%L2NORM0 = dnrm2(this%NEQ, this%D, 1)
    !
    ! -- CHECK FOR EXACT SOLUTION
    itmax = this%ITER1
    IF (this%L2NORM0 == DZERO) THEN
      itmax = 0
      ICNVG = 1
    END IF
    !
    ! -- SOLUTION BY THE CONJUGATE GRADIENT METHOD
    IF (this%ILINMETH == 1) THEN
      CALL ims_base_cg(ICNVG, itmax, innerit, &
                       this%NEQ, this%NJA, this%NIAPC, this%NJAPC, &
                       this%IPC, this%ICNVGOPT, this%NORTH, &
                       this%DVCLOSE, this%RCLOSE, this%L2NORM0, &
                       this%EPFACT, this%IA0, this%JA0, this%A0, &
                       this%IAPC, this%JAPC, this%APC, &
                       this%X, this%RHS, this%D, this%P, this%Q, this%Z, &
                       this%NJLU, this%IW, this%JLU, &
                       NCONV, CONVNMOD, CONVMODSTART, &
                       CACCEL, summary)
      !
      ! -- SOLUTION BY THE BICONJUGATE GRADIENT STABILIZED METHOD
    ELSE IF (this%ILINMETH == 2) THEN
      CALL ims_base_bcgs(ICNVG, itmax, innerit, &
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
                         CACCEL, summary)
    END IF
    !
    ! -- BACK PERMUTE AMAT, SOLUTION, AND RHS
    IF (this%IORD /= 0) THEN
      CALL dperm(this%NEQ, this%A0, this%JA0, this%IA0, &
                 this%AMAT, this%JA, this%IA, &
                 this%IORDER, this%ID, 1)
      CALL dvperm(this%NEQ, this%X, this%IORDER)
      CALL dvperm(this%NEQ, this%RHS, this%IORDER)
    END IF
    !
    ! -- UNSCALE PROBLEM
    IF (this%ISCL .NE. 0) THEN
      CALL ims_base_scale(1, this%ISCL, &
                          this%NEQ, this%NJA, this%IA, this%JA, &
                          this%AMAT, this%X, this%RHS, &
                          this%DSCALE, this%DSCALE2)
    END IF
    !
    ! -- SET IMS INNER ITERATION NUMBER (IN_ITER) TO NUMBER OF
    !       IMSLINEAR INNER ITERATIONS (innerit)
    IN_ITER = innerit
  end SUBROUTINE imslinear_ap

END MODULE IMSLinearModule
