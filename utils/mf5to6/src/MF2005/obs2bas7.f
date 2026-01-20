      MODULE OBSBASMODULE
         INTEGER, SAVE,POINTER ::ITS,NH,MAXM,MOBS,IUHOBSV,IDRY,JDRY
         INTEGER, SAVE,POINTER ::IPRT
         double precision,    SAVE,POINTER ::HOBDRY
         INTEGER, SAVE, DIMENSION(:,:), POINTER ::NDER
         INTEGER, SAVE, DIMENSION(:,:), POINTER ::MLAY
         INTEGER, SAVE, DIMENSION(:),   POINTER ::IOFF
         INTEGER, SAVE, DIMENSION(:),   POINTER ::JOFF
         INTEGER, SAVE, DIMENSION(:),   POINTER ::IHOBWET
         REAL,    SAVE, DIMENSION(:),   POINTER ::H
         REAL,    SAVE, DIMENSION(:),   POINTER ::HOBS
         REAL,    SAVE, DIMENSION(:),   POINTER ::TOFF
         double precision,    SAVE, DIMENSION(:),   POINTER ::ROFF
         double precision,    SAVE, DIMENSION(:),   POINTER ::COFF
         REAL,    SAVE, DIMENSION(:),   POINTER ::OTIME
         REAL,    SAVE, DIMENSION(:,:), POINTER ::PR
         REAL,    SAVE, DIMENSION(:,:), POINTER ::RINT
         CHARACTER*12,SAVE,DIMENSION(:),POINTER ::OBSNAM
         integer, save, dimension(:),   pointer :: irefspd
         integer, save, dimension(:),   pointer :: nlayer
       TYPE OBSBASTYPE
         INTEGER,  POINTER    ::ITS,NH,MAXM,MOBS,IUHOBSV,IDRY,JDRY
         INTEGER,  POINTER    ::IPRT
         double precision,     POINTER    ::HOBDRY
         INTEGER,  DIMENSION(:,:), POINTER ::NDER
         INTEGER,  DIMENSION(:,:), POINTER ::MLAY
         INTEGER,  DIMENSION(:),   POINTER ::IOFF
         INTEGER,  DIMENSION(:),   POINTER ::JOFF
         INTEGER,  DIMENSION(:),   POINTER ::IHOBWET
         REAL,     DIMENSION(:),   POINTER ::H
         REAL,     DIMENSION(:),   POINTER ::HOBS
         REAL,     DIMENSION(:),   POINTER ::TOFF
         double precision,     DIMENSION(:),   POINTER ::ROFF
         double precision,     DIMENSION(:),   POINTER ::COFF
         REAL,     DIMENSION(:),   POINTER ::OTIME
         REAL,     DIMENSION(:,:), POINTER ::PR
         REAL,     DIMENSION(:,:), POINTER ::RINT
         CHARACTER*12,DIMENSION(:),POINTER ::OBSNAM
         integer,  dimension(:),   pointer :: irefspd
         integer,  dimension(:),   pointer :: nlayer
      END TYPE
      TYPE(OBSBASTYPE), SAVE :: OBSBASDAT(10)
      END MODULE OBSBASMODULE
      
C  NDER(1,n) -- Observation layer
C  NDER(2,n) -- Observation row
C  NDER(3,n) -- Observation column
C  NDER(4,n) -- Observation time step
C  NDER(5,n) -- Observation number for computing observation as a head change
C  MLAY(MAXM,MOBS) -- Layer numbers for multilayer observations
C  IOFF(NH) -- Row offset for neighboring cell for interpolation
C  JOFF(NH) -- Column offset for neighboring cell for interpolation
C  IHOBWET(NH) -- Flag for observation -- 1 for wet, and -1 for dry
C  H(NH) -- Simulated value
C  HOBS(NH) -- Observed value
C  TOFF(NH) -- Fractional offset between time steps
C  ROFF(NH) -- Fractional offset from center of cell in Y direction (between rows)
C  COFF(NH) -- Fractional offset from center of cell in X direction (between columns)
C  OTIME(NH) -- Observation time in model time units
C  PR(MAXM,MOBS) -- Fractional value for each layer of multilayer observations
C  RINT(4,NH) -- Interpolation coefficients for the 4 nodes surrounding observation
C  OBSNAM(NH) -- Observation name

      SUBROUTINE OBS2BAS7AR(IUHDOB,IGRID,needPreproc)
C     ******************************************************************
C     INITIALIZE AND READ VARIABLES FOR HEAD OBSERVATIONS
C     ******************************************************************
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY: NCOL,NROW,NLAY,DELR,DELC,
     1                  NPER,NSTP,PERLEN,TSMULT,ISSFLG,IOUT,ITRSS
      USE OBSBASMODULE
      use SimPHMFModule, only: ustop
      use utl7module, only: URDCOM, URWORD
C
      logical, intent(inout) :: needPreproc
      CHARACTER*200 LINE
      double precision :: dum
C     ------------------------------------------------------------------
C
      needPreproc = .false.
C1------ALLOCATE AND INITIALIZE TIME STEP COUNTER FOR USE BY ANY
C1------OBSERVATION PACKAGE.
      ALLOCATE(ITS)
      ITS=0
      IF(IUHDOB.LE.0) GO TO 700
C
C2------ALLOCATE OTHER SCALARS IF HEAD OBSERVATIONS ARE BEING SPECIFIED.
      ALLOCATE(NH,MAXM,MOBS,IUHOBSV,IDRY,JDRY,IPRT)
      ALLOCATE(HOBDRY)
C
C3------DEFINE CONSTANTS
      IDRY = 0
      JDRY = 0
      ZERO=0.0
      ONEN=-1.0
      ML=0
      IERR=0
C
C4------WRITE VERSION.
      WRITE (IOUT,14) IUHDOB
   14 FORMAT (/,' OBS2BAS7 -- HEAD OBSERVATIONS, ',
     &        'VERSION 2.0, 2/28/2006',/,' INPUT READ FROM UNIT ',I3)
C
C5------READ & PRINT ITEM 1 OF THE HOB INPUT FILE
      CALL URDCOM(IUHDOB,IOUT,LINE)
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NH,DUM,IOUT,IUHDOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MOBS,DUM,IOUT,IUHDOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MAXM,DUM,IOUT,IUHDOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUHOBSV,DUM,IOUT,IUHDOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,HOBDRY,IOUT,IUHDOB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,IDUM,DUM,IOUT,IUHDOB)
      IPRT=1
      IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
        IPRT=0
        WRITE(IOUT,*) 'NOPRINT option for HEAD OBSERVATIONS'
      END IF
      IF (MAXM.EQ.1) THEN
        WRITE (IOUT,17)
   17   FORMAT (/,' MAXM CAN NOT EQUAL 1 -- STOP EXECUTION')
        CALL USTOP(' ')
      ENDIF
      WRITE (IOUT,19) NH, MOBS, MAXM
   19 FORMAT (/,
     &     ' NUMBER OF HEADS....................................:',I5,/,
     &     '   NUMBER OF MULTILAYER HEADS.......................:',I5,/,
     &     '   MAXIMUM NUMBER OF LAYERS FOR MULTILAYER HEADS....:',I5)
      IF(NH.LE.0) THEN
         WRITE(IOUT,*) ' NH LESS THAN OR EQUAL TO 0'
         CALL USTOP(' ')
      END IF
      IF(IUHOBSV.GT.0) THEN
         WRITE(IOUT,21) IUHOBSV
   21    FORMAT(1X,
     1      'HEAD OBSERVATIONS WILL BE SAVED ON UNIT............:',I5)
      ELSE
         WRITE(IOUT,22)
   22    FORMAT(1X,'HEAD OBSERVATIONS WILL NOT BE SAVED IN A FILE')
      END IF
      WRITE(IOUT,23) HOBDRY
   23 FORMAT(1X,'SIMULATED EQUIVALENT HEAD AT DRY CELLS WILL BE:',
     1            1P,1G15.6)
C
C6------ALLOCATE ARRAY DATA.
      ALLOCATE(NDER(5,NH))
      ALLOCATE(IOFF(NH))
      ALLOCATE(JOFF(NH))
      ALLOCATE(IHOBWET(NH))
      ALLOCATE(OBSNAM(NH))
      ALLOCATE(H(NH))
      ALLOCATE(HOBS(NH))
      ALLOCATE(TOFF(NH))
      ALLOCATE(OTIME(NH))
      ALLOCATE(ROFF(NH))
      ALLOCATE(COFF(NH))
      roff = 0.0d0
      coff = 0.0d0
      ALLOCATE(RINT(4,NH))
      IF(MOBS.GT.0 .AND. MAXM.GT.1) THEN
        ALLOCATE(MLAY(MAXM,MOBS))
        ALLOCATE(PR(MAXM,MOBS))
      ELSE
        ALLOCATE(MLAY(1,1))
        ALLOCATE(PR(1,1))
      END IF
      allocate(irefspd(nh))
      allocate(nlayer(nh))
C
C7------INITIALIZE OTIME, SIMULATED EQUIVALENT HEAD, NDER(5,n), and IHOBWET.
      DO N = 1, NH
        OTIME(N) = ONEN
        H(N) = ZERO
        NDER(5,N)=0
        IHOBWET(N)=1
        irefspd(n) = 0
        nlayer(n) = 0
      enddo
C
C8------READ ITEM 2
      READ (IUHDOB,*) TOMULTH
C
C9------WRITE ITEM 2 AND TITLE FOR OBSERVATION TIMES.
      IF(IPRT.NE.0) WRITE (IOUT,530) TOMULTH
  530 FORMAT (/,' OBSERVED HEAD DATA -- TIME OFFSETS ARE',
     &' MULTIPLIED BY: ',G12.5,//,
     &20X,'REFER.',/,
     &7X,'OBSERVATION',2X,'STRESS',4X,'TIME',/,
     &2X,'OBS#    NAME',6X,'PERIOD',3X,'OFFSET    OBSERVATION')
C
C10-----INITIALIZE N, WHICH IS THE COUNT OF THE NUMBER OF OBSERVATIONS
C10-----AS THEY ARE READ.
      N=0
C
C11-----READ NAME, LOCATION, TIME, AND OBSERVED VALUE (ITEM 3)
   60 N=N+1
      READ (IUHDOB,*) OBSNAM(N), (NDER(I,N),I=1,3), IREFSP, TOFFSET,
     &                 ROFF(N), COFF(N), HOBS(N)
      if (roff(n) /= 0.0d0 .or. coff(n) /= 0.0d0) then
        needPreproc = .true.
      endif
      irefspd(n) = irefsp
      IF(IPRT.NE.0) WRITE (IOUT,535) N,OBSNAM(N),IREFSP,TOFFSET,HOBS(N)
  535 FORMAT (1X,I5,1X,A12,2X,I4,2X,G11.4,1X,G11.4)
C
C11A----FOR SINGLE-TIME OBSERVATION (IREFSP>0), CALL UOBSTI TO DETERMINE
C11A----WHEN OBSERVATION OCCURS.
      IF (IREFSP.GE.0) THEN
        CALL UOBSTI(OBSNAM(N),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,
     &                NDER(4,N),PERLEN,TOFF(N),TOFFSET,TOMULTH,TSMULT,
     &                0,OTIME(N))
      END IF
C
C12-----CHECK ROW AND COLUMN LOCATION.
      I = NDER(2,N)
      J = NDER(3,N)
      IF (J.LE.0 .OR. J.GT.NCOL .OR. I.LE.0 .OR. I.GT.NROW) THEN
        WRITE (IOUT,550) N
  550   FORMAT (' FOR OBS',I5,' ROW OR COLUMN NUMBER INVALID -- ',
     &        'STOP EXECUTION (OBS2BAS7HRP)',/)
        IERR = 1
      ENDIF
C
C13-----Check if multi-layer
      IF(NDER(1,N).GE.0) THEN
C
C13A----SINGLE LAYER -- CHECK FOR VALID LAYER.
        IF (NDER(1,N).LE.0 .OR. NDER(1,N).GT.NLAY) THEN
          WRITE (IOUT,265) NDER(1,N)
  265     FORMAT (' FOR OBS ',I5,' LAYER INVALID -- STOP EXECUTION',
     &        ' (OBS2BAS7AR)',/)
          IERR = 1
        END IF
        nlayer(n) = 1
      ELSE
C
C13B----MULTI-LAYER -- CHECK LIMITS AND READ THE LAYERS AND PROPORTIONS.
         needPreproc = .true.
         NL=-NDER(1,N)
         nlayer(n) = nl
         ML = ML + 1
         IF(ML.GT.MOBS) THEN
           WRITE (IOUT,565)
  565 FORMAT (/,' NUMBER OF MULTILAYER OBSERVATIONS EXCEEDS MOBS -- ',
     &        'STOP EXECUTION (OBS2BAS7AR)',/)
           CALL USTOP(' ')
         END IF
         IF(NL.GT.MAXM) THEN
           WRITE(IOUT,620) NL
  620      FORMAT(/,1X,'ERROR: VALUE ENTERED FOR MAXM IN HOB FILE IS',
     &    ' SMALLER THAN THE MAXIMUM NUMBER',/,' OF LAYERS',
     &    ' IN A MULTILAYER HEAD OBSERVATION, WHICH IS ',
     &    I3,' -- INCREASE MAXM.')
           CALL USTOP(' ')
         END IF
         DO M=1,MAXM
           MLAY(M,ML) = 0
         enddo
         READ (IUHDOB,*) (MLAY(M,ML),PR(M,ML),M=1,NL)
         IF(IPRT.NE.0) WRITE(IOUT,540) (MLAY(M,ML),PR(M,ML),M=1,NL)
  540  FORMAT (5X,'MULTIPLE LAYERS AND PROPORTIONS :',5(I5,',',F5.2,3X))
C
C
C13C----CHECK LAYER NUMBERS AND ADD PROPORTIONS FOR MULTILAYER
C13C----OBSERVATION WELLS.
        TPR=ZERO
        DO K = 1,NL
          KK = MLAY(K,ML)
          TPR = TPR + PR(K,ML)
          IF (KK.LE.0 .OR. KK.GT.NLAY) THEN
            WRITE (IOUT,265) N
            IERR = 1
          ENDIF
        enddo
C
C13D----CHECK SUM OF PROPORTIONS FOR MULTILAYER OBS WELLS
        IF (ABS(1.-TPR).GT..02) THEN
          WRITE (IOUT,560) N
  560 FORMAT (/,' FOR OBS',I5,' MULTILAYER PROPORTIONS DO NOT SUM ',
     &        'TO 1.0 -- STOP EXECUTION (OBS2BAS7AR)',/)
          IERR = 1
        ENDIF
      END IF
C
C14-----CALCULATE INTERPOLATION COEFFICIENTS FOR THE LOCATION.
      CALL SOBS2BAS7HIA(N,ML)
C
C15-----CHECK FOR MULTI-TIME
      NT=-IREFSP
      IF(NT.GT.0) THEN
C
C15A----READ FLAG FOR USING TEMPORAL CHANGES IN HEAD (ITEM 5)
        READ (IUHDOB,*) ITT
        IF(IPRT.NE.0) WRITE (IOUT,515) ITT
  515   FORMAT (2X,'TRANSIENT DATA AT THIS LOCATION, ITT =',I4)
        IF (ITT.NE.1 .AND. ITT.NE.2) THEN
          WRITE (IOUT,575) N
  575     FORMAT (' FOR OBS',I5,
     &     ' ITT MUST = 1 OR 2 -- STOP EXECUTION (OBS2BAS7HRP)',/)
          CALL USTOP(' ')
        ENDIF
C
C15B----LOOP THROUGH THE TIMES
        NBASE=N
        DO J=1,NT
        IF(J.NE.1) THEN
C
C15B1---DUPLICATE THE LOCATION INFORMATION FOR THE OBSERVATIONS AT THE
C15B1---SAME LOCATION.
          N=N+1
          nlayer(n) = nlayer(n-1)
          IF(N.GT.NH) THEN
            WRITE(IOUT,127)
  127       FORMAT(1X,/,1X,'ABORTING BECAUSE THERE ARE MORE HEAD',
     1                       ' OBSERVATIONS THAN SPECIFIED BY NH')
            CALL USTOP(' ')
          END IF
          DO I = 1, 3
            NDER(I,N) = NDER(I,N-1)
          enddo
          ROFF(N) = ROFF(N-1)
          COFF(N) = COFF(N-1)
          IOFF(N) = IOFF(N-1)
          JOFF(N) = JOFF(N-1)
          DO I = 1, 4
            RINT(I,N) = RINT(I,N-1)
          enddo
          IF (NDER(1,N-1).LT.0) THEN
            ML1 = ML
            ML = ML + 1
            DO M = 1, MAXM
              PR(M,ML) = PR(M,ML1)
              MLAY(M,ML) = MLAY(M,ML1)
            enddo
          ENDIF
        END IF 
C
C15B2---READ ONE MULTI-TIME OBSERVATION.
        READ (IUHDOB,*) OBSNAM(N), IREFSP, TOFFSET, HOBS(N)
        irefspd(n) = -irefsp
        IF(ITT.EQ.2 .AND. J.NE.1) THEN
           HOBS(N)=HOBS(N)-HOBS(NBASE)
           NDER(5,N)=NBASE
        END IF
C
C15B3---WRITE ONE MULTI-TIME OBSERVATION.
        IF(IPRT.NE.0) WRITE (IOUT,535) N,OBSNAM(N),IREFSP,TOFFSET,
     1                                 HOBS(N)
        CALL UOBSTI(OBSNAM(N),IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,
     &              NDER(4,N),PERLEN,TOFF(N),TOFFSET,TOMULTH,
     &              TSMULT,0,OTIME(N))
        IF(J.EQ.NT .AND. IPRT.NE.0) WRITE (IOUT,570)
  570   FORMAT (' ')
        enddo
      END IF
C
C16-----READ ANOTHER OBSERVATION (ITEM 3) IF THERE ARE STILL MORE OBSERVATIONS.
      IF(N.LT.NH) GO TO 60
C
C17-----DONE READING HEAD OBSERVATIONS.
C17-----PRINT TABLE SHOWING LOCATION OF OBSERVATIONS.
      IF(IPRT.NE.0) THEN
        WRITE (IOUT,590)
  590   FORMAT (/,53X,'HEAD CHANGE',/,54X,'REFERENCE',/,
     &8X,'OBSERVATION',19X,'ROW',5X,'COL    OBSERVATION',/,
     &2X,'OBS#',5X,'NAME',7X,'LAY  ROW  COL  OFFSET  OFFSET',3X,
     &'(IF > 0)')
        DO N = 1, NH
          WRITE (IOUT,600) N, OBSNAM(N), (NDER(I,N),I=1,3), ROFF(N),
     &                   COFF(N), NDER(5,N)
  600     FORMAT (1X,I5,2X,A12,2X,I3,2(1X,I4),2(2X,F6.3),3X,I6)
        enddo
      END IF
C
C18-----IF ERROR OCCURRED ABOVE, PRINT MESSAGE AND STOP.
      IF (IERR.GT.0) THEN
        WRITE(IOUT,610)
  610 FORMAT(/,1X,'ERROR: SEARCH ABOVE FOR ERROR MESSAGE(S)',/,1X,
     &'STOP EXECUTION -- (OBS2BAS7AR)')
        CALL USTOP(' ')
      ENDIF
C
C19-----RETURN.
  700 CALL SOBS2BAS7PSV(IUHDOB,IGRID)
      RETURN
      END SUBROUTINE OBS2BAS7AR
      
      SUBROUTINE SOBS2BAS7HIA(N,ML)
C     ******************************************************************
C     CALCULATE INTERPOLATION COEFFICIENTS FOR LOCATING OBSERVED HEADS
C     ASSUMING ALL CELLS ARE ACTIVE.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL, ONLY: NCOL,NROW,NLAY,DELR,DELC,IBOUND,HNEW,STRT,
     1                  NPER,NSTP,PERLEN,TSMULT,ISSFLG,IOUT
      USE OBSBASMODULE
C     ------------------------------------------------------------------
C
      K = NDER(1,N)
      IF (K.LT.0) K = MLAY(1,ML)
      I = NDER(2,N)
      J = NDER(3,N)
      I1 = I + 1
      J1 = J + 1
      IOFF(N) = 1
      JOFF(N) = 1
      IF (ROFF(N).LT.0.) THEN
        I1 = I - 1
        IOFF(N) = -1
      ENDIF
      IF (COFF(N).LT.0.) THEN
        J1 = J - 1
        JOFF(N) = -1
      ENDIF
      IF (I1.GE.1 .AND. I1.LE.NROW) IBI = 1
      IF (J1.GE.1 .AND. J1.LE.NCOL) IBJ = 1
      IF (I1.GE.1 .AND. I1.LE.NROW .AND. J1.GE.1 .AND. J1.LE.NCOL)
     &    IBIJ = 1
      IF (I1.LT.1 .OR. I1.GT.NROW) THEN
        ROFF(N) = 0.
        IBI = 0
      ENDIF
      IF (J1.LT.1 .OR. J1.GT.NCOL) THEN
        COFF(N) = 0.
        IBJ = 0
      ENDIF
      IF (I1.LT.1 .OR. I1.GT.NROW .OR. J1.LT.1 .OR. J1.GT.NCOL) IBIJ = 0
C
      CALL SOBS2BAS7HBF(COFF(N),DELC,DELR,I,I1,IBI,IBIJ,IBJ,IOFF(N),
     1                  J,J1,JOFF(N),NCOL,NROW,RINT(:,N),ROFF(N))
C
      RETURN
      END SUBROUTINE SOBS2BAS7HIA
      
      SUBROUTINE SOBS2BAS7HIB(NDER,COFF,ROFF,DELR,DELC,IBOUND,NCOL,NROW,
     &                      NLAY,RINT,JOFF,IOFF,MLAY)
C     ******************************************************************
C     CALCULATE INTERPOLATION COEFFICIENTS FOR LOCATING OBSERVED HEADS
C     USING CURRENT IBOUND VALUES.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION NDER(5), IBOUND(NCOL,NROW,NLAY),
     &          RINT(4)
      double precision :: coff, roff, delr(ncol), delc(nrow)
C     ------------------------------------------------------------------
C
      K = NDER(1)
      IF (K.LT.0) K = MLAY
      I = NDER(2)
      J = NDER(3)
      I1 = I + 1
      J1 = J + 1
      IOFF = 1
      JOFF = 1
      IF (ROFF.LT.0.) THEN
        I1 = I - 1
        IOFF = -1
      ENDIF
      IF (COFF.LT.0.) THEN
        J1 = J - 1
        JOFF = -1
      ENDIF
      IF (I1.GE.1 .AND. I1.LE.NROW) IBI = IBOUND(J,I1,K)
      IF (J1.GE.1 .AND. J1.LE.NCOL) IBJ = IBOUND(J1,I,K)
      IF (I1.GE.1 .AND. I1.LE.NROW .AND. J1.GE.1 .AND. J1.LE.NCOL)
     &    IBIJ = IBOUND(J1,I1,K)
      IF (I1.LT.1 .OR. I1.GT.NROW) THEN
        ROFF = 0.
        IBI = 0
      ENDIF
      IF (J1.LT.1 .OR. J1.GT.NCOL) THEN
        COFF = 0.
        IBJ = 0
      ENDIF
      IF (I1.LT.1 .OR. I1.GT.NROW .OR. J1.LT.1 .OR. J1.GT.NCOL) IBIJ = 0
C
      CALL SOBS2BAS7HBF(COFF,DELC,DELR,I,I1,IBI,IBIJ,IBJ,IOFF,J,J1,JOFF,
     &                NCOL,NROW,RINT,ROFF)
C
      RETURN
      END SUBROUTINE SOBS2BAS7HIB
      
      SUBROUTINE SOBS2BAS7HBF(COFF,DELC,DELR,I,I1,IBI,IBIJ,IBJ,IOFF,J,
     &                      J1,JOFF,NCOL,NROW,RINT,ROFF)
C     ******************************************************************
C     CALCULATE BASIS FUNCTIONS FOR INTERPOLATING OBSERVED HEADS.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      double precision :: coff, roff, delr, delc
      DIMENSION DELR(NCOL), DELC(NROW), RINT(4)
C     ------------------------------------------------------------------
      A=0.
C
C1------MOVE OBSERVATION TO NODE IF CLOSE TO NODE OR IF NEIGHBORS ARE
C1------NO FLOW
      IF ((ABS(ROFF).LT..001.AND.ABS(COFF).LT..001) .OR.
     &    (ABS(ROFF).LT..001.AND.IBJ.EQ.0) .OR.
     &    (ABS(COFF).LT..001.AND.IBI.EQ.0) .OR. (IBI.EQ.0.AND.IBJ.EQ.0))
     &    THEN
        IOFF = 0
        JOFF = 0
        DO IR = 1, 4
          RINT(IR) = .25
        enddo
        RETURN
      ENDIF
C
C2------CALCULATE CONSTANTS
      IF (ABS(ROFF).GE..001) THEN
        DC = (DELC(I)+DELC(I1))/2.
        DCF = ABS(ROFF)*DELC(I)
      ENDIF
      IF (ABS(COFF).GE..001) THEN
        DR = (DELR(J)+DELR(J1))/2.
        DRF = ABS(COFF)*DELR(J)
      ENDIF
      IF (ABS(ROFF).GE..001 .AND. ABS(COFF).GE..001) A = 1/(DC*DR)
C
C3------LINEAR INTERPOLATION
      IF (ABS(ROFF).LT..001 .OR. (IBI.EQ.0.AND.IBIJ.EQ.0)) THEN
        IOFF = 0
        RINT(1) = 0.5*(1.-DRF/DR)
        RINT(2) = 0.5*DRF/DR
        RINT(3) = RINT(1)
        RINT(4) = RINT(2)
C
      ELSEIF (ABS(COFF).LT..001 .OR. (IBJ.EQ.0.AND.IBIJ.EQ.0)) THEN
        JOFF = 0
        RINT(1) = 0.5*(1.-DCF/DC)
        RINT(2) = RINT(1)
        RINT(3) = 0.5*DCF/DC
        RINT(4) = RINT(3)
C
C4------CALCULATE BASIS FUNCTIONS FOR INTERPOLATION ON A RECTANGLE
      ELSEIF (IBJ.NE.0 .AND. IBI.NE.0 .AND. IBIJ.NE.0) THEN
        RINT(3) = A*(DR-DRF)*DCF
        RINT(4) = A*DRF*DCF
        RINT(2) = A*DRF*(DC-DCF)
        RINT(1) = A*(DR-DRF)*(DC-DCF)
C
C5------CALCULATE BASIS FUNCTIONS FOR INTERPOLATION ON A TRIANGLE
      ELSEIF (IBJ.EQ.0) THEN
        RINT(1) = A*(DR*DC-DR*DCF)
        RINT(2) = 0.0
        RINT(3) = A*(DR*DCF-DC*DRF)
        RINT(4) = A*(DC*DRF)
C
      ELSEIF (IBI.EQ.0) THEN
        RINT(1) = A*(DR*DC-DC*DRF)
        RINT(4) = A*(DR*DCF)
        RINT(2) = A*(DC*DRF-DR*DCF)
        RINT(3) = 0.0
C
      ELSEIF (IBIJ.EQ.0) THEN
        RINT(1) = A*(DR*DC-DC*DRF-DR*DCF)
        RINT(3) = A*(DR*DCF)
        RINT(2) = A*(DC*DRF)
        RINT(4) = 0.0
      ENDIF
C
C6------
      RETURN
      END SUBROUTINE SOBS2BAS7HBF
      
      SUBROUTINE UOBSTI(ID,IOUT,ISSFLG,ITRSS,NPER,NSTP,IREFSP,NUMTS,
     &                  PERLEN,TOFF1,TOFFSET,TOMULT,TSMULT,ITR1ST,
     &                  OBSTIME)
      use SimPHMFModule, only: ustop
C     ******************************************************************
C     ASSIGN OBSERVATION TIME STEP (NUMTS) AND TOFF GIVEN REFERENCE
C     STRESS PERIOD (IREFSP), OBSERVATION-TIME OFFSET (TOFFSET), AND
C     TIME-OFFSET MULTIPLIER (TOMULT)
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*(*) ID
      INTEGER IREFSP, ISSFLG, ITR1ST, NPER, NSTP, NUMTS
      REAL DELT, ENDTIME, TIME, TOFF1, TOFFMULT, TOFFSET,
     &     TOMULT
      double precision :: perlen, TSMULT
      DIMENSION NSTP(NPER), PERLEN(NPER), TSMULT(NPER), ISSFLG(NPER)
C     ------------------------------------------------------------------
      ZERO=0.
C
C1------ENSURE THAT SPECIFIED REFERENCE STRESS PERIOD IS VALID
      IF (IREFSP.LT.1 .OR. IREFSP.GT.NPER) THEN
        WRITE(IOUT,505) IREFSP
 505    FORMAT(/,' REFERENCE STRESS PERIOD (IREFSP) WAS SPECIFIED AS ',
     &  I5,', BUT IT MUST BE',/,
     &  ' BETWEEN 1 AND NPER (OF THE DISCRETIZATION INPUT FILE)',/,
     &  ' -- STOP EXECUTION (UOBSTI)')
        CALL USTOP(' ')
      ENDIF
C
C2------ENSURE THAT TOFFSET IS NOT NEGATIVE
      IF (TOFFSET.LT.ZERO) THEN
        WRITE(IOUT,510) TRIM(ID)
 510    FORMAT(/,' TOFFSET IS NEGATIVE FOR OBSERVATION "',A,
     &  '" -- STOP EXECUTION (UOBSTI)')
        CALL USTOP(' ')
      ENDIF
C
C3------FIND NUMBER OF TIME STEPS PRECEDING REFERENCE STRESS PERIOD
      NUMTS = 0
      OBSTIME = ZERO
      IF (IREFSP.GT.1) THEN
        DO I = 1, IREFSP-1
          NUMTS = NUMTS + NSTP(I)
          OBSTIME = OBSTIME + PERLEN(I)
        enddo
      ENDIF
C
C Note that variables TIME and ENDTIME are relative to the reference stress
C  period. Variable OBSTIME is relative to the start of the simulation.
      TIME = ZERO
C
C4------USE TOMULT TO CONVERT TOFFSET TO MODEL-TIME UNITS (ASSUMES THAT
C4------USER HAS DEFINED TOMULT CORRECTLY)
      TOFFMULT = TOFFSET*TOMULT
      !
      ! Catch special case where, in MF2005, for a steady-state stress period,
      ! TOMULTH*TOFFSET < end of stress period indicates observation is for
      ! the steady solution. In MF6, the time must be specified for the end of
      ! the steady-state stress period.
      if (ISSFLG(irefsp) == 1) then
        ! Stress period IREFSP is steady-state
        if (toffmult < perlen(irefsp)) then
          toffmult = perlen(irefsp)
        endif
      endif
      !
      OBSTIME = OBSTIME + TOFFMULT
C
C5------FIND STRESS PERIOD IN WHICH OBSERVATION TIME FALLS.
C5------LOOP THROUGH STRESS PERIODS STARTING AT REFERENCE STRESS PERIOD.
C5------TOFF1 IS OBSERVATION TIME IN TIME STEP, AS A FRACTION OF THE TIME
C5------STEP. NUMTS IS THE NUMBER OF THE TIME STEP PRECEDING THE TIME STEP
C5------IN WHICH THE OBSERVATION TIME OCCURS.
      DO I = IREFSP, NPER
        ENDTIME = TIME+PERLEN(I)
        IF (ENDTIME.GE.TOFFMULT) THEN
C
C6------FIND TIME STEP PRECEDING OBSERVATION TIME
C         CALCULATE LENGTH OF FIRST TIME STEP IN CURRENT STRESS PERIOD
          DELT = PERLEN(I)/FLOAT(NSTP(I))
          IF (TSMULT(I).NE.1.) DELT = PERLEN(I)*(1.-TSMULT(I))/
     &                                (1.-TSMULT(I)**NSTP(I))
C
C7------LOOP THROUGH TIME STEPS
          DO J = 1, NSTP(I)
            ENDTIME = TIME+DELT
            IF (ENDTIME.GE.TOFFMULT) THEN
              IF (ISSFLG(I).NE.0 .OR. ITRSS.EQ.0) THEN
C
C8------STEADY-STATE TIME STEP.
C8------SET NUMTS AS THE START OF THE NEXT TIME STEP.  TELL USER UNLESS
C8------THE OBSERVATION TIME IS THE END OF THE TIME STEP
                IF(TOFFMULT.LT.ENDTIME) THEN
                  WRITE(IOUT,33)
   33   FORMAT(1X,'Observation within a steady-state time step has',
     1       ' been moved to the end of the time step.')
                END IF
                TOFF1 = 1.0
              ELSE
C
C9------Transient time step.
C9A-----CALCULATE TOFF1 AS FRACTION OF TIME-STEP DURATION
                TOFF1 = (TOFFMULT-TIME)/DELT
C
C9B-----CHECK FOR INITIAL TRANSIENT TIME STEP
                IF (NUMTS.EQ.0) THEN
                  IF(ITR1ST.EQ.1) THEN
C
C9C-----ITR1ST IS A FLAG THAT INDICATES IF, FOR A CERTAIN
C9C-----OBSERVATION TYPE, THE OBSERVATION TIME CAN BE BEFORE THE
C9C-----END OF AN INITIAL TRANSIENT TIME STEP.  ITR1ST=0
C9C-----INDICATES THAT THE OBS. TIME CAN BE IN INITIAL TRANSIENT
C9C-----TIME STEP.  ITR1ST=1 INDICATES IT CAN'T.
                    WRITE(IOUT,37)
   37      FORMAT(1X,'The observation is in the first time step of the',
     1  ' simulation, but the observation type does not allow this.',/
     2 1X,'The observation is being moved to the end of the time step.')
                    TOFF1 = 1.0
                  ELSE
C
C9D-----STOP IF THE OBSERVATION IS AT THE BEGINNING OF AN INITIAL TRANSIENT
C9D-----TIME STEP.
                    IF(TOFFMULT.EQ.ZERO) THEN
                      WRITE(IOUT,38)
   38     FORMAT(1X,'An observation cannot be placed at the very',
     1 ' beginning of the simulation if the first period is transient.')
                      CALL USTOP(' ')
                    END IF
                  END IF
                ENDIF
              ENDIF
              GOTO 80
            ENDIF
            TIME = TIME+DELT
            DELT = DELT*TSMULT(I)
            NUMTS = NUMTS+1
          enddo
        ELSE
          NUMTS = NUMTS+NSTP(I)
          TIME = TIME+PERLEN(I)
        ENDIF
      enddo
C
C10-----ALLOW FOR ROUND-OFF ERROR, SO THAT OBSERVATION TIMES SPECIFIED
C10-----AT THE EXACT END OF THE SIMULATION ARE NOT FLAGGED AS ERRORS
      TOLERANCE = 1.0E-6*PERLEN(NPER)
      TDIFF = TOFFMULT-TIME
      IF (TDIFF.LT.TOLERANCE) THEN
        TOFF1 = 1.0
        NUMTS=NUMTS-1
      ELSE
        WRITE(IOUT,500) ID
 500    FORMAT(/,' TIME SPECIFIED FOR OBSERVATION "',A,
     & '" IS AFTER END OF SIMULATION',/,' -- STOP EXECUTION (UOBSTI)')
        CALL USTOP(' ')
      ENDIF
C
C11-----The Time step and interpolation coefficient have been determined.
 80   CONTINUE
      RETURN
      END SUBROUTINE UOBSTI
      
      SUBROUTINE UOBSSV(IUOBSSV,NOBS,H,HOBS,OBSNAM,LABEL)
C     ******************************************************************
C     SAVE OBSERVATIONS TO A DISK FILE
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION H(NOBS),HOBS(NOBS)
      CHARACTER*(*) OBSNAM(NOBS)
C     ------------------------------------------------------------------
C
      IF(IUOBSSV.GT.0) THEN
C
C1------WRITE LABEL IF "LABEL" IS NOT 0
        IF(LABEL.NE.0) WRITE(IUOBSSV,18)
   18   FORMAT('"SIMULATED EQUIVALENT"',3X,'"OBSERVED VALUE"',
     1       4X,'"OBSERVATION NAME"')
C
C2------WRITE OBSERVATIONS
        DO N=1,NOBS
          WRITE(IUOBSSV,28) H(N),HOBS(N),OBSNAM(N)
   28     FORMAT(1X,1P,E19.11,E20.11,2X,A)
        enddo
      END IF
C
C3------RETURN
      RETURN
      END SUBROUTINE UOBSSV
      
      SUBROUTINE OBS2BAS7DA(IUHDOB,IGRID)
C  Deallocate OBSBAS memory
      USE OBSBASMODULE
C
      CALL SOBS2BAS7PNT(IUHDOB,IGRID)
      DEALLOCATE(ITS)
      IF(IUHDOB.LE.0) RETURN
C
      DEALLOCATE(NH)
      DEALLOCATE(MAXM)
      DEALLOCATE(MOBS)
      DEALLOCATE(IUHOBSV)
      DEALLOCATE(IDRY)
      DEALLOCATE(JDRY)
      DEALLOCATE(IPRT)
      DEALLOCATE(HOBDRY)
      DEALLOCATE(NDER)
      DEALLOCATE(MLAY)
      DEALLOCATE(IOFF)
      DEALLOCATE(JOFF)
      DEALLOCATE(IHOBWET)
      DEALLOCATE(H)
      DEALLOCATE(HOBS)
      DEALLOCATE(TOFF)
      DEALLOCATE(ROFF)
      DEALLOCATE(COFF)
      DEALLOCATE(OTIME)
      DEALLOCATE(PR)
      DEALLOCATE(RINT)
      DEALLOCATE(OBSNAM)
      deallocate(irefspd)
      deallocate(nlayer)
C
      RETURN
      END SUBROUTINE OBS2BAS7DA
      
      SUBROUTINE SOBS2BAS7PNT(IUHDOB,IGRID)
C  Change OBSBAS data to a different grid.
      USE OBSBASMODULE
C
      ITS=>OBSBASDAT(IGRID)%ITS
      IF(IUHDOB.LE.0) RETURN
C
      NH=>OBSBASDAT(IGRID)%NH
      MAXM=>OBSBASDAT(IGRID)%MAXM
      MOBS=>OBSBASDAT(IGRID)%MOBS
      IUHOBSV=>OBSBASDAT(IGRID)%IUHOBSV
      IDRY=>OBSBASDAT(IGRID)%IDRY
      JDRY=>OBSBASDAT(IGRID)%JDRY
      IPRT=>OBSBASDAT(IGRID)%IPRT
      HOBDRY=>OBSBASDAT(IGRID)%HOBDRY
      NDER=>OBSBASDAT(IGRID)%NDER
      MLAY=>OBSBASDAT(IGRID)%MLAY
      IOFF=>OBSBASDAT(IGRID)%IOFF
      JOFF=>OBSBASDAT(IGRID)%JOFF
      IHOBWET=>OBSBASDAT(IGRID)%IHOBWET
      H=>OBSBASDAT(IGRID)%H
      HOBS=>OBSBASDAT(IGRID)%HOBS
      TOFF=>OBSBASDAT(IGRID)%TOFF
      ROFF=>OBSBASDAT(IGRID)%ROFF
      COFF=>OBSBASDAT(IGRID)%COFF
      OTIME=>OBSBASDAT(IGRID)%OTIME
      PR=>OBSBASDAT(IGRID)%PR
      RINT=>OBSBASDAT(IGRID)%RINT
      OBSNAM=>OBSBASDAT(IGRID)%OBSNAM
      irefspd=>OBSBASDAT(igrid)%irefspd
      nlayer=>OBSBASDAT(igrid)%nlayer
C
      RETURN
      END SUBROUTINE SOBS2BAS7PNT
      
      SUBROUTINE SOBS2BAS7PSV(IUHDOB,IGRID)
C  Save OBSBAS data for a grid.
      USE OBSBASMODULE
C
C
      OBSBASDAT(IGRID)%ITS=>ITS
      IF(IUHDOB.LE.0) RETURN
C
      OBSBASDAT(IGRID)%NH=>NH
      OBSBASDAT(IGRID)%MAXM=>MAXM
      OBSBASDAT(IGRID)%MOBS=>MOBS
      OBSBASDAT(IGRID)%IUHOBSV=>IUHOBSV
      OBSBASDAT(IGRID)%IDRY=>IDRY
      OBSBASDAT(IGRID)%JDRY=>JDRY
      OBSBASDAT(IGRID)%IPRT=>IPRT
      OBSBASDAT(IGRID)%HOBDRY=>HOBDRY
      OBSBASDAT(IGRID)%NDER=>NDER
      OBSBASDAT(IGRID)%MLAY=>MLAY
      OBSBASDAT(IGRID)%IOFF=>IOFF
      OBSBASDAT(IGRID)%JOFF=>JOFF
      OBSBASDAT(IGRID)%IHOBWET=>IHOBWET
      OBSBASDAT(IGRID)%H=>H
      OBSBASDAT(IGRID)%HOBS=>HOBS
      OBSBASDAT(IGRID)%TOFF=>TOFF
      OBSBASDAT(IGRID)%ROFF=>ROFF
      OBSBASDAT(IGRID)%COFF=>COFF
      OBSBASDAT(IGRID)%OTIME=>OTIME
      OBSBASDAT(IGRID)%PR=>PR
      OBSBASDAT(IGRID)%RINT=>RINT
      OBSBASDAT(IGRID)%OBSNAM=>OBSNAM
      OBSBASDAT(igrid)%irefspd=>irefspd
      OBSBASDAT(igrid)%nlayer=>nlayer
C
      RETURN
      END SUBROUTINE SOBS2BAS7PSV
