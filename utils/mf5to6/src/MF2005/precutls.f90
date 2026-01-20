!
!  This Module must be compiled with default real KIND = 4 bytes 
!                    and default double precision KIND = 8 bytes
!

MODULE PRECUTLSMOD

  use InputOutputModule, only: same_word
  use SimPHMFModule, only: store_error, store_note, store_warning, ustop
  use UtilitiesModule, only: findcell

  PRIVATE
  PUBLIC :: GET_BINARY_HEAD_DATASET, GET_BUDGET_3D_ARRAY, &
            GET_BINARY_FILE_LENGTH, BUDGETPRECISION, &
            HEADPRECISION
  PUBLIC :: AUXDIM
  INTEGER, PARAMETER :: AUXDIM = 20
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IBUF
  REAL, ALLOCATABLE, DIMENSION(:) :: BUF
  REAL, ALLOCATABLE, DIMENSION(:,:) :: BUFF2D
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BUFF3D
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DBUF
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DBUF3D
  REAL, DIMENSION(6) :: QV
  DOUBLE PRECISION, DIMENSION(6) :: QVD
  
CONTAINS
  
  INTEGER FUNCTION BUDGETPRECISION(IU)
    ! Determine single or double precision file type for a MODFLOW
    ! budget file:  0=unrecognized, 1=single, 2=double.
    IMPLICIT NONE
    ! Argument-list variables
    INTEGER, INTENT(IN) :: IU
    ! Local variables
    INTEGER :: ICELL, ICODE, IERR, IPREC, KPER, KSTP, &
               N, NC, NCOL, NL, NLAY, NODES, NLST, NR, NROW
    REAL :: DELT, PERTIM, TOTIM, VAL
    DOUBLE PRECISION :: DELTD, PERTIMD, TOTIMD, VALD
    CHARACTER(LEN=16) :: TEXT1, TEXT2
    LOGICAL :: MOREDATA
    !
    !  Default is unrecognized file
    IPREC=0
    !
    !  SINGLE check
    READ(IU,ERR=100,END=100) KSTP,KPER,TEXT1,NCOL,NROW,NLAY
    ICODE=0
    IF (NLAY.LT.0) THEN
      NLAY=-NLAY
      READ(IU,ERR=50,END=50) ICODE,DELT,PERTIM,TOTIM
    ENDIF
    IF (NCOL.LT.1 .OR. NROW.LT.1 .OR. NLAY.LT.1) GO TO 100
    IF (NCOL.GT.100000000 .OR. NROW.GT.100000000 .OR. &
                        NLAY.GT.100000000) GO TO 100
    IF (NCOL*NROW.GT.100000000 .OR. NCOL*NLAY.GT.100000000 .OR. &
                    NROW*NLAY.GT.100000000) GO TO 100
    ALLOCATE(BUFF3D(NCOL,NROW,NLAY))
    ALLOCATE(DBUF3D(NCOL,NROW,NLAY))
    NODES=NCOL*NROW*NLAY
    !
    !  Read data depending on ICODE.  
    IF (ICODE.EQ.0 .OR. ICODE.EQ.1) THEN
      READ(IU,ERR=50,END=50) BUFF3D
    ELSEIF (ICODE.EQ.2) THEN
      READ(IU,ERR=50,END=50) NLST
      IF (NLST.LT.0) GO TO 50
      IF (NLST.GT.0) THEN
        DO N=1,NLST
          READ(IU,END=50,ERR=50) ICELL,VAL
          IF (ICELL.LE.0 .OR. ICELL.GT.NODES) GO TO 50
        ENDDO
      ENDIF
    ELSEIF (ICODE==3) THEN
      CALL READ_2D_ARRAY_LAYER(IU,1,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 50
    ELSEIF (ICODE==4) THEN
      CALL READ_2D_ARRAY_L1(IU,1,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 50
    ELSEIF (ICODE==5) THEN
      CALL READ_LIST_AUX(IU,1,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 50
    ELSE
      GO TO 100
    ENDIF
    !
    !  Read 2nd header and check for valid type.
    MOREDATA = .FALSE.
    READ(IU,ERR=50,END=40) KSTP,KPER,TEXT2
    MOREDATA = .TRUE.
    IF (VALID_BUDGET_TEXT(TEXT1) .AND. VALID_BUDGET_TEXT(TEXT2)) THEN
      IPREC = 1
    ENDIF
    40 CONTINUE
    IF (.NOT. MOREDATA) IPREC = 1 ! EOF indicates file contains a single, valid data set
    IF (IPREC==1) GOTO 100
    !
    !  DOUBLE check
    50 REWIND(IU)
    READ(IU,ERR=100,END=100) KSTP,KPER,TEXT1,NC,NR,NL
    ICODE=0
    IF (NL.LT.0) THEN
      NL=-NL
      READ(IU,ERR=100,END=100) ICODE,DELTD,PERTIMD,TOTIMD
    ENDIF
    !
    !  Read data depending on ICODE.
    IF (ICODE.EQ.0 .OR. ICODE.EQ.1) THEN
      READ(IU,ERR=100,END=100) DBUF3D
    ELSEIF (ICODE.EQ.2) THEN
      READ(IU,ERR=100,END=100) NLST
      IF (NLST.LT.0) GO TO 100
      IF (NLST.GT.0) THEN
        DO N=1,NLST
          READ(IU,END=100,ERR=100) ICELL,VALD
          IF (ICELL.LE.0 .OR. ICELL.GT.NODES) GO TO 100
        ENDDO
      ENDIF
    ELSEIF (ICODE==3) THEN
      CALL READ_2D_ARRAY_LAYER(IU,2,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 100
    ELSEIF (ICODE==4) THEN
      CALL READ_2D_ARRAY_L1(IU,2,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 100
    ELSEIF (ICODE==5) THEN
      CALL READ_LIST_AUX(IU,2,NCOL,NROW,NLAY,IERR,DBUF3D)
      IF (IERR/=0) GOTO 100
    ELSE
      GO TO 100
    ENDIF
    !
    !  Read 2nd header and check for valid type.
    MOREDATA = .FALSE.
    READ(IU,ERR=100,END=90) KSTP,KPER,TEXT2
    MOREDATA = .TRUE.
    IF (VALID_BUDGET_TEXT(TEXT1) .AND. VALID_BUDGET_TEXT(TEXT2)) THEN
      IPREC = 2
    ENDIF
    90 CONTINUE
    IF (.NOT. MOREDATA) IPREC = 2
    !
    100 REWIND(IU)
    BUDGETPRECISION = IPREC
    IF (ALLOCATED(BUFF3D)) DEALLOCATE(BUFF3D)
    IF (ALLOCATED(DBUF3D)) DEALLOCATE(DBUF3D)
    RETURN
  END FUNCTION BUDGETPRECISION
  !-------------------------------------------------------------------
  INTEGER FUNCTION HEADPRECISION(IU)
    ! Return 1 if head or drawdown file opened on unit IU is single
    ! precision.  Return 2 if it's double precision.  Return value
    ! of -1 means precision cannot be determined.
    IMPLICIT NONE
    ! Argument
    INTEGER, INTENT(IN) :: IU
    ! Local variables
    INTEGER :: I, ILAY, ISTAT, J, KP, KS, NCOL, NCOL2, NROW, NROW2
    CHARACTER(LEN=16) :: TEXT
    REAL :: PERTIMS, TOTIMS
    REAL, ALLOCATABLE, DIMENSION(:,:) :: VALS
    DOUBLE PRECISION :: PERTIMD, TOTIMD
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VALD
    !
    HEADPRECISION = -1 ! Unknown
    if (iu > 0) then
      REWIND(IU)
      !
      ! First assume single precision
      READ(IU,ERR=20,END=20)KS,KP,PERTIMS,TOTIMS,TEXT,NCOL,NROW,ILAY
      IF (NCOL>0 .AND. NROW>0) THEN
        IF (NCOL*NROW<100000000) THEN
          ALLOCATE(VALS(NCOL,NROW),STAT=ISTAT)
          IF (ISTAT==0) THEN
            DO I=1,NROW
              READ(IU,ERR=20,END=20)(VALS(J,I),J=1,NCOL)
            ENDDO
            READ(IU,ERR=20,END=10)KS,KP,PERTIMS,TOTIMS,TEXT,NCOL2,NROW2,ILAY
            IF (NCOL2==NCOL .AND. NROW2==NROW) THEN
              HEADPRECISION = 1
            ELSE 
              GOTO 20
            ENDIF
          ELSE
            GOTO 20
          ENDIF
        ELSE
          GOTO 20
        ENDIF
      ELSE 
        GOTO 20
      ENDIF
      !
      10 CONTINUE
      ! If binary head file contains only 1 2-D head array, READ of 2nd header 
      ! record will fail with EOF error.  Since other READs must have been 
      ! successful, file is single precision
      HEADPRECISION = 1
      !
      ! Error reading as single precision
      20 CONTINUE
      !
      IF (HEADPRECISION<0) THEN
        ! Retry, assuming double precision
        REWIND(IU)
        READ(IU,ERR=40,END=40)KS,KP,PERTIMD,TOTIMD,TEXT,NCOL,NROW,ILAY
        IF (NCOL>0 .AND. NROW>0) THEN
          IF (NCOL*NROW<100000000) THEN
            ALLOCATE(VALD(NCOL,NROW),STAT=ISTAT)
            IF (ISTAT==0) THEN
              DO I=1,NROW
                READ(IU,ERR=40,END=40)(VALD(J,I),J=1,NCOL)
              ENDDO
              READ(IU,ERR=40,END=30)KS,KP,PERTIMD,TOTIMD,TEXT,NCOL2,NROW2,ILAY
              IF (NCOL2==NCOL .AND. NROW2==NROW) THEN
                HEADPRECISION = 2
              ELSE
                GOTO 40
              ENDIF
            ELSE
              GOTO 40
            ENDIF
          ELSE
            GOTO 40
          ENDIF
        ELSE
          GOTO 40
        ENDIF
        30 CONTINUE
        ! If binary head file contains only 1 2-D head array, READ of 2nd header 
        ! record will fail.  Since other READs must have been successful, file
        ! is double precision
        HEADPRECISION = 2
      ENDIF      
      !
      ! Deallocate arrays and exit
      40 CONTINUE
      IF (ALLOCATED(VALS)) DEALLOCATE(VALS)
      IF (ALLOCATED(VALD)) DEALLOCATE(VALD)
      REWIND(IU)
    else
      ! Return flag indicating IU <= 0
      headprecision = -2
    endif
    !
    RETURN
  END FUNCTION HEADPRECISION
  !-------------------------------------------------------------------
  SUBROUTINE GET_BINARY_HEAD_DATASET(IU,PRECISION,KSTP,KPER,TEXT, &
                                      NCOL,NROW,ILAY,DHEADS,OK)
    ! Read binary head data written by ULASAV
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, NCOL, NROW, PRECISION
    INTEGER, INTENT(OUT) :: KSTP, KPER, ILAY
    CHARACTER(LEN=16), INTENT(OUT) :: TEXT
    DOUBLE PRECISION, DIMENSION(NCOL,NROW), INTENT(OUT) :: DHEADS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: I, J, NC, NR
    REAL :: PERTIMS, TOTIMS
    DOUBLE PRECISION :: PERTIMD, TOTIMD
    LOGICAL :: OKLOCAL
    !
    OKLOCAL = .TRUE.
    DHEADS = 0.0D0
    SELECT CASE (PRECISION)
    CASE (1)
      IF (.NOT. ALLOCATED(BUFF2D)) ALLOCATE(BUFF2D(NCOL,NROW))
      READ(IU,ERR=900,END=950)KSTP,KPER,PERTIMS,TOTIMS,TEXT,NC,NR,ILAY
      DO I=1,NR
        READ(IU,ERR=900,END=950)(BUFF2D(J,I),J=1,NC)
        DO J=1,NC
          DHEADS(J,I) = BUFF2D(J,I)
        ENDDO
      ENDDO
    CASE (2)
      READ(IU,ERR=900,END=950)KSTP,KPER,PERTIMD,TOTIMD,TEXT,NC,NR,ILAY
      DO I=1,NR
        READ(IU,ERR=900,END=950)(DHEADS(J,I),J=1,NC)
      ENDDO
    END SELECT
    TEXT = ADJUSTL(TEXT)
    !
    ! Normal return
    RETURN
    !
    ! Error handling
    900 CONTINUE
    CALL store_error('Error extracting simulated head value(s)')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
    950 CONTINUE
    CALL store_error('Binary head file missing simulated value(s)')
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE GET_BINARY_HEAD_DATASET
  !-------------------------------------------------------------------
  SUBROUTINE GET_BUDGET_3D_ARRAY(IU,IPREC,NCOL,NROW,NLAY, &
                                 QUANTITY,NCB,KPER,KSTP,DBUF3D,OK)
    ! Read a data set from binary cell-by-cell budget file written by 
    ! (1) UBUDSV (3D array), or (2) UBDSV1 (3D array with time data), 
    ! or (3) UBDSV2 and UBDSVA (list), or (4) UBDSV3 (2D array with 
    ! optional 2D array of layer numbers), or (5) UBDSV4 and UBDSVB 
    ! (list).  Return data in 3D array DBUF3D.  If IPREC==1, the 
    ! budget file is assumed to be single precision; if IPREC==2, it 
    ! is assumed to be double precision.
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NCOL, NLAY, NROW
    CHARACTER(LEN=*), INTENT(OUT) :: QUANTITY ! TEXT associated with quantity of interest
    INTEGER, INTENT(OUT) :: KPER, KSTP
    DOUBLE PRECISION, DIMENSION(NCOL,NROW,NLAY), INTENT(OUT) :: DBUF3D
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: IERR, IST, NC, NCB, NL, NR
    REAL :: DELTLOCAL, PERTIM, TOTIM
    DOUBLE PRECISION :: DELTDLOCAL, PERTIMD, TOTIMD
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=16) :: TEXT
    !
    OKLOCAL = .TRUE.
    DBUF3D = 0.0D0
    !
    ! Read a header record from unformatted budget file
    READ(IU,ERR=900,END=950,IOSTAT=IST) KSTP,KPER,TEXT,NC,NR,NL
    QUANTITY = ADJUSTL(TEXT)
    IF (IST.NE.0) THEN
      OKLOCAL = .FALSE.
      GOTO 900
    ENDIF
    ! Negative NLAY indicates COMPACT BUDGET
    IF (NL.LT.0) THEN
      NL=-NL
      ! Read the extra record that defines the compact type
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) NCB,DELTLOCAL,PERTIM,TOTIM
      CASE (2)
        READ(IU,END=950,ERR=900) NCB,DELTDLOCAL,PERTIMD,TOTIMD
      END SELECT
    ELSE
      NCB = 1
    ENDIF
    !
    SELECT CASE (NCB)
    CASE (1)
      ! Full 3-D array (many packages when budget file is not COMPACT)
      CALL READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (2)
      !  List without auxiliary values (SFR)
      CALL READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (3)
      !  2-D array plus a layer indicator array (ETS, EVT, RCH)
      CALL READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (4)
      !  2-D array for layer 1 (ETS, EVT, RCH)
      CALL READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    CASE (5)
      !  List with auxiliary values (WEL, DRN, RIV, GHB)
      CALL READ_LIST_AUX(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    END SELECT
    IF (IERR==1) THEN
      GOTO 900
    ELSEIF (IERR==2) THEN
      GOTO 950
    ENDIF
    !
    ! Normal return
    RETURN
    !
    ! Error handling
    900 CONTINUE
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    CALL store_error('Error extracting data from cell-by-cell budget file')
    RETURN
    950 CONTINUE
    OKLOCAL = .FALSE.
    IF (OK) OK = OKLOCAL
    CALL store_error('Binary cell-by-cell budget file missing simulated value(s)')
    RETURN
  END SUBROUTINE GET_BUDGET_3D_ARRAY
  !-------------------------------------------------------------------
  SUBROUTINE READ_3D_ARRAY(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, J, K, N, NRCL
    !
    IERR = 0
    NRCL = NC*NR*NL
    SELECT CASE (IPREC)
    CASE (1)
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NC*NR*NL))
      READ(IU,END=950,ERR=900) (BUF(I),I=1,NRCL)
    CASE (2)
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NC*NR*NL))
      READ(IU,END=950,ERR=900) (DBUF(I),I=1,NRCL)
    CASE DEFAULT
      GOTO 900
    END SELECT
    IF (.TRUE.) THEN
      ! Copy values from BUF to 3-D array
      N = 0
      DO K=1,NL
        DO I=1,NR
          DO J=1,NC
            N = N+1
            SELECT CASE (IPREC)
            CASE (1)
              DBUF3D(J,I,K) = REAL(BUF(N),8)
            CASE (2)
              DBUF3D(J,I,K) = DBUF(N)
            END SELECT
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_3D_ARRAY
  !-------------------------------------------------------------------
  SUBROUTINE READ_LIST(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, ICRL, J, K, L, NLIST
    REAL :: Q
    DOUBLE PRECISION :: QD
    !
    IERR = 0
    READ(IU,END=950,ERR=900) NLIST
    DO L=1,NLIST
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) ICRL,Q
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + REAL(Q,8)
      CASE (2)
        READ(IU,END=950,ERR=900) ICRL,QD
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + QD
      END SELECT
    ENDDO
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_LIST
  !-------------------------------------------------------------------
  SUBROUTINE READ_LIST_SFR(IU,IPREC,NSTRM,QD,IERR)
    ! Read NLIST (written by UBDSV2) and reach-by-reach data
    ! written by SFR using UBDSVA
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NSTRM
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NSTRM), INTENT(OUT) :: QD
    ! Local variables
    INTEGER :: ICRL, L, NLIST
    REAL :: Q
    !
    QD = 0.0D0
    IERR = 0
    READ(IU,END=950,ERR=900) NLIST
    DO L=1,NLIST
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) ICRL,Q
        QD(L) = REAL(Q,8)
      CASE (2)
        READ(IU,END=950,ERR=900) ICRL,QD(L)
      END SELECT
    ENDDO
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_LIST_SFR
  !-------------------------------------------------------------------
  SUBROUTINE READ_LIST_AUX(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, ICRL, J, K, L, N, NLIST, NVAL
    CHARACTER(LEN=16) :: TEXT
    !
    IERR = 0
    READ(IU,END=950,ERR=900) NVAL
    IF (NVAL.GT.1) THEN
      READ(IU,END=950,ERR=900) (TEXT,N=2,NVAL)
    ENDIF
    READ(IU,END=950,ERR=900) NLIST
    DO L=1,NLIST
      SELECT CASE (IPREC)
      CASE (1)
        READ(IU,END=950,ERR=900) ICRL,(QV(N),N=1,NVAL)
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + REAL(QV(1),8)
      CASE (2)
        READ(IU,END=950,ERR=900) ICRL,(QVD(N),N=1,NVAL)
        CALL FINDCELL(ICRL,NR,NC,NL,I,J,K)
        IF (K.EQ.0) GOTO 900
        DBUF3D(J,I,K) = DBUF3D(J,I,K) + QVD(1)
      END SELECT
    ENDDO
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_LIST_AUX
  !-------------------------------------------------------------------
  SUBROUTINE READ_2D_ARRAY_LAYER(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, J, K, N, NRC
    !
    IERR = 0
    NRC = NC*NR
    IF (.NOT. ALLOCATED(IBUF)) ALLOCATE(IBUF(NRC))
    READ(IU,END=950,ERR=900) (IBUF(N),N=1,NRC)
    SELECT CASE (IPREC)
    CASE (1)
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NRC))
      READ(IU,END=950,ERR=900) (BUF(N),N=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,IBUF(N)) = REAL(BUF(N),8)
      ENDDO
    CASE (2)
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NRC))
      READ(IU,END=950,ERR=900) (DBUF(I),I=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,IBUF(N)) = DBUF(N)
        ENDDO
    END SELECT
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_2D_ARRAY_LAYER
  !-------------------------------------------------------------------
  SUBROUTINE READ_2D_ARRAY_L1(IU,IPREC,NC,NR,NL,IERR,DBUF3D)
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, NC, NR, NL
    INTEGER, INTENT(OUT) :: IERR
    DOUBLE PRECISION, DIMENSION(NC,NR,NL), INTENT(INOUT) :: DBUF3D
    ! Local variables
    INTEGER :: I, J, K, N, NRC
    !
    IERR = 0
    NRC = NC*NR
    SELECT CASE (IPREC)
    CASE (1)
      IF (.NOT. ALLOCATED(BUF)) ALLOCATE(BUF(NRC))
      READ(IU,END=950,ERR=900) (BUF(I),I=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,1) = REAL(BUF(N),8)
      ENDDO
    CASE (2)
      IF (.NOT. ALLOCATED(DBUF)) ALLOCATE(DBUF(NRC))
      READ(IU,END=950,ERR=900) (DBUF(I),I=1,NRC)
      DO N=1,NRC
        CALL FINDCELL(N,NR,NC,1,I,J,K)
        IF (I.EQ.0) GOTO 900        
        DBUF3D(J,I,1) = DBUF(N)
      ENDDO
    END SELECT
    !
    ! Normal return    
    RETURN
    !
    ! Error handling
    900 CONTINUE
    IERR = 1  ! Unidentified error
    RETURN
    !
    950 CONTINUE
    IERR = 2  ! End of file
    RETURN
  END SUBROUTINE READ_2D_ARRAY_L1
  !-------------------------------------------------------------------
  SUBROUTINE READ_NCB_5(IU,IPREC,MAXDIM,AUXDIM,AUXTXT,AUXVALS, &
                        MAXAUXVALUE,AUXVAR,FNAME,AUXINDX,NAUX, &
                        NLIST,ICRLS,QVALS,OK)
    ! Invoke when NCB = 5 (data written by UBDSV4 and UBDSVB).
    ! Read list of values with or without auxiliary values, starting
    ! with NAUX+1 (written by UBDSV4), followed (optionally) by
    ! AUXTXT (names of AUX variables), followed by NLIST and NLIST 
    ! repetitions of ICRL, Q, and (optionally) AUX values (written
    ! by UBDSVB).  Names of AUX variables are returned in AUXTXT, 
    ! ICRL values are returned in ICRLS, and Q values are returned 
    ! in QVALS.  NAUX read from binary file is returned.
    IMPLICIT NONE
    ! Arguments
    INTEGER, INTENT(IN) :: IU, IPREC, MAXDIM, AUXDIM
    CHARACTER(LEN=16) :: AUXTXT(AUXDIM)
    DOUBLE PRECISION, DIMENSION(AUXDIM,MAXDIM) :: AUXVALS
    ! MAXAUXVALUE is the maximum aux-value read from mmproc.in.
    INTEGER, INTENT(IN) :: MAXAUXVALUE
    ! AUXVAR IS name of AUX variable that identifies boundary 
    ! (e.g. drain) of interest.
    CHARACTER(LEN=*), INTENT(IN) :: AUXVAR, FNAME
    ! AUXINDX is the index in AUXVALS that
    ! corresponds to AUX variable AUXVAR.
    INTEGER, INTENT(OUT) :: AUXINDX
    INTEGER, INTENT(OUT) :: NAUX, NLIST
    INTEGER, DIMENSION(MAXDIM), INTENT(INOUT) :: ICRLS
    DOUBLE PRECISION, DIMENSION(MAXDIM),   INTENT(INOUT) :: QVALS
    LOGICAL, INTENT(INOUT) :: OK
    ! Local variables
    INTEGER :: II, N, NAUXP1
    REAL :: AUXVAL, Q
    LOGICAL :: OKLOCAL
    CHARACTER(LEN=400) :: ERRMSG
    !
    OKLOCAL = .TRUE.
    AUXINDX = 0
    !
    IF (OKLOCAL) THEN
      ! Read remaining data written by UBDSV4
      READ(IU) NAUXP1
      NAUX = NAUXP1 - 1
!      READ(IU) NAUX
      IF (NAUX>0) READ(IU) (AUXTXT(N),N=1,NAUX)
      IF (MAXAUXVALUE>-1) THEN
        ! Assign AUXINDX
        IF (NAUX>0) THEN
          DO N=1,NAUX
            IF (same_word(AUXTXT(N),AUXVAR)) THEN
              AUXINDX = N
              EXIT
            ENDIF
          ENDDO
        ENDIF
        IF (AUXINDX==0) THEN
          ERRMSG = 'Required AUX variable ("' // TRIM(AUXVAR) // &
              '") not found in binary file "' // TRIM(FNAME) // &
              '". Are you missing "COMPACT BUDGET AUXILIARY"' // &
              ' in the Output Control input file?'
          CALL store_error(ERRMSG)
          OKLOCAL = .FALSE.
          GOTO 900
        ENDIF
      ENDIF
      READ(IU) NLIST
      !
      IF (NLIST>MAXDIM) THEN
        ERRMSG = 'MAXDIM dimension too small in READ_NCB_5'
        CALL store_error(ERRMSG)
        OKLOCAL = .FALSE.
        GOTO 900
      ENDIF
      !
      ! Read data written by UBDSVB
      IF (IPREC==1) THEN
        DO II=1,NLIST
          READ(IU) ICRLS(II)
          READ(IU) Q
          QVALS(II) = REAL(Q,8)
          IF(NAUX.GT.0) THEN
             DO N=1,NAUX
               READ(IU)AUXVAL
               AUXVALS(N,II) = REAL(AUXVAL,8)
             ENDDO
          ELSE
          ENDIF
        ENDDO
      ELSEIF (IPREC==2) THEN
        DO II=1,NLIST
          READ(IU) ICRLS(II),QVALS(II)
          IF(NAUX.GT.0) THEN
            READ(IU) (AUXVALS(N,II),N=1,NAUX)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
    900 CONTINUE
    IF (OK) OK = OKLOCAL
    RETURN
  END SUBROUTINE READ_NCB_5
  !-------------------------------------------------------------------
  INTEGER FUNCTION GET_BINARY_FILE_LENGTH(IU)
    ! Return length of a binary file, in bytes.
    ! Assume unit IU has been opened (return -1 if not).
    ! Rewind IU before return.
    IMPLICIT NONE
    ! Argument
    INTEGER, INTENT(IN) :: IU
    ! Local variables
    INTEGER :: K
    CHARACTER(LEN=1) :: BYT
    LOGICAL :: LOP, OK
    !
    INQUIRE(UNIT=IU,OPENED=LOP)
    IF (LOP) THEN
      REWIND(IU)
      K = 0
      OK = .TRUE.
      DO WHILE (OK)
        READ(IU,END=10)BYT
        K = K + 1
      ENDDO
      10 CONTINUE
      REWIND(IU)
      GET_BINARY_FILE_LENGTH = K
    ELSE
      GET_BINARY_FILE_LENGTH = -1
    ENDIF
    !
    RETURN
  END FUNCTION GET_BINARY_FILE_LENGTH
  !-------------------------------------------------------------------
  LOGICAL FUNCTION VALID_BUDGET_TEXT(TEXT)
    ! Determine if TEXT is a valid text entry in a cell-by-cell budget
    ! file.  Try to make it support all MODFLOW packages.
    IMPLICIT NONE
    ! Argument
    CHARACTER(LEN=*), INTENT(IN) :: TEXT
    ! Local variables
    CHARACTER(LEN=16) :: LOCALTEXT, LEFTTEXT
    LOGICAL :: VALID
    INTEGER :: I, TEXTSDIM, TEXTLEN
    PARAMETER (TEXTSDIM=30)
    CHARACTER(LEN=16), DIMENSION(TEXTSDIM) :: TEXTS
    DATA TEXTS/'CONSTANT HEAD   ', &
               'FLOW RIGHT FACE ', &
               'FLOW FRONT FACE ', &
               'FLOW LOWER FACE ', &
               'STORAGE         ', &
               'DRAINS          ', &
               'DRAINS (DRT)    ', &
               'ET SEGMENTS     ', &
               'ET              ', &
               'SPECIFIED FLOWS ', &
               'HEAD DEP BOUNDS ', &
               'INTERBED STORAGE', &
               'LAKE SEEPAGE    ', &
               'MNW             ', &
               'MNW2            ', &
               'RECHARGE        ', &
               'RESERV. LEAKAGE ', &
               'RIVER LEAKAGE   ', &
               'STREAM LEAKAGE  ', &
               'STREAMFLOW OUT  ', &
               'STREAM LEAKAGE  ', &
               'STREAM FLOW OUT ', &
               'INST. IB STORAGE', &
               'DELAY IB STORAGE', &
               'SWR LEAKAGE     ', &
               'SWR GWET        ', &
               'UZF RECHARGE    ', &
               'GW ET           ', &
               'SURFACE LEAKAGE ', &
               'WELLS           '/               
    !
    ! If TEXT is blank, it's invalid
    IF (TEXT==' ') THEN
      VALID_BUDGET_TEXT = .FALSE.
      RETURN
    ENDIF
    !
    VALID = .FALSE.
    TEXTLEN = LEN(TEXT)
    !
    LOCALTEXT = TEXT
    ! Left-justify text
    LEFTTEXT = ADJUSTL(LOCALTEXT)
    !
    ! Iteratively compare left-justified text with valid values
    COMPARE: DO I=1,TEXTSDIM
      IF (same_word(LEFTTEXT,TEXTS(I))) THEN
        VALID = .TRUE.
        EXIT COMPARE
      ENDIF
    ENDDO COMPARE
    !
    VALID_BUDGET_TEXT = VALID
    RETURN
  END FUNCTION VALID_BUDGET_TEXT
  
END MODULE PRECUTLSMOD
