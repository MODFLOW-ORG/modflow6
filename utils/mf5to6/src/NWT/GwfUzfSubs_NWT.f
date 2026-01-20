      module GwfUzfSubs

      use InputOutputModule, only: URWORD
      use ArrayReadersMF5Module, only: ReadArray
      use SimPHMFModule, only: ustop
      use utl7module, only: URDCOM !, UBUDSV, UBDSVB, UBDSV3
        
      contains
C
C-------SUBROUTINE GWF2UZF1AR
      SUBROUTINE GWF2UZF1AR(In, Iunitbcf, Iunitlpf,  
     +                      Iunitupw, Igrid)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR UNSATURATED FLOW, RECHARGE, AND ET
C     READ AND CHECK VARIABLES THAT REMAIN CONSTANT
!--------REVISED FOR MODFLOW-2005 RELEASE 1.9, FEBRUARY 6, 2012
!rgn------REVISION NUMBER CHANGED TO BE CONSISTENT WITH NWT RELEASE
!rgn------NEW VERSION NUMBER 1.0.9:  July 1, 2014
C     ******************************************************************
      USE GWFUZFMODULE
      USE GLOBAL,       ONLY: NCOL, NROW, NLAY, IOUT, ITRSS, ISSFLG, 
     +                        DELR, DELC, IBOUND, LBOTM, BOTM
      USE GLOBAL,       ONLY: ITMUNI, LENUNI, LAYHDT
      USE GWFLPFMODULE, ONLY: SCLPF=>SC2
      USE GWFBCFMODULE, ONLY: SC1, SC2, LAYCON
*      USE GWFHUFMODULE, ONLY: SC2HUF
      USE GWFUPWMODULE, ONLY: SC2UPW
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ARGUMENTS
C     ------------------------------------------------------------------
      INTEGER In, Iunitbcf, Iunitlpf, Iunitupw, Igrid
C     ------------------------------------------------------------------
C     LOCAL VARIABLES
C     ------------------------------------------------------------------
      DOUBLE PRECISION test
      INTEGER istart, istop, lloc, ivol, numactive, ic, ir
      INTEGER ibndflg, ichld, iflgbnd, igage, igunit, irhld, isyflg, 
     +        iuzcol, iuzflg, iuzlay, iuzopt, iuzrow, l, ncck, ncth, 
     +        nlth, nrck, nrnc, nrth, i, icheck, kkrch, k, NPP, MXVL,
     +        llocsave
      double precision :: r, sy, fkmin, fkmax, range, finc, thick
      CHARACTER(LEN=200) line
      CHARACTER(LEN=24) aname(8)
      DATA aname(1)/' AREAL EXTENT OF UZ FLOW'/
      DATA aname(2)/' ROUTING OVERLAND RUNOFF'/
      DATA aname(3)/' SATURATED WATER CONTENT'/
      DATA aname(4)/'   INITIAL WATER CONTENT'/
      DATA aname(5)/'    BROOKS-COREY EPSILON'/
      DATA aname(6)/'    SATURATED VERTICAL K'/
      DATA aname(7)/'    UZ CELL BOTTOM ELEV.'/
      DATA aname(8)/'  RESIDUAL WATER CONTENT'/
C     ------------------------------------------------------------------
      Version_uzf =
     +'$Id: gwf2uzf1_NWT.f 4071 2014-07-01 23:30:24Z rniswon $'
      ALLOCATE(NUMCELLS, TOTCELLS, Iseepsupress, IPRCNT)
      Iseepsupress = 0   ! Iseepsupress =1 means seepout not calculated
      NUMCELLS = NCOL*NROW
      TOTCELLS = NUMCELLS*NLAY
      IPRCNT = 0
      ALLOCATE (LAYNUM(NCOL,NROW))
      ALLOCATE (NUZTOP, IUZFOPT, IRUNFLG, IETFLG, IUZM)
      ALLOCATE (IUZFCB1, IUZFCB2, NTRAIL, NWAV, NSETS, IUZFB22, IUZFB11)
      ALLOCATE (NUZGAG, NUZGAGAR, NUZCL, NUZRW, TOTRUNOFF)
      ALLOCATE (SURFDEP,IGSFLOW, RTSOLUTE)
      ALLOCATE (ITHTIFLG, ITHTRFLG, IETBUD)
      allocate (NeedUzfWaterMover)
      ITHTIFLG = 0
      ITHTRFLG = 0
      IETBUD = 0
C
C1------IDENTIFY PACKAGE AND INITIALIZE.
      WRITE (IOUT, 9001) In
 9001 FORMAT (1X, /' UZF1 -- UNSATURATED FLOW PACKAGE, VERSION 1.0.5', 
     +        ', 10/01/2011', /, 9X, 'INPUT READ FROM UNIT', I3)
      CALL URDCOM(In, IOUT, line)
! Check for alternate input.
      CALL UPARLSTAL(IN,IOUT,LINE,NPP,MXVL)
      lloc = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SPECIFYTHTR') THEN
        ITHTRFLG = 1
        WRITE(iout,*)
        WRITE(IOUT,'(A)')' RESIDUAL WATER CONTENT (THTR) WILL BE READ ',
     +                ' AND USED FOR THE FIRST TRANSIENT STRESS PERIOD'
        WRITE(iout,*)
        llocsave = lloc
      ELSE
        WRITE(iout,*)
        WRITE(IOUT,'(A)')' RESIDUAL WATER CONTENT (THTR) WILL BE ',
     +              ' CALCULATED AS THE DIFFERENCE BETWEEN THTS AND SY'
        WRITE(iout,*)
        llocsave = 1
      END IF
      lloc = llocsave
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SPECIFYTHTI') THEN
         ITHTIFLG = 1
         WRITE(iout,*)
         WRITE(IOUT,'(A)')' INITIAL WATER CONTENT (THTI) WILL BE READ ',
     +                ' FOR THE FIRST SS OR TR STRESS PERIOD'
         WRITE(iout,*)
         llocsave = lloc
      ELSE
         WRITE(iout,*)
         WRITE(IOUT,'(A)') ' INITIAL WATER CONTENT (THTI) WILL BE ',
     +                 ' CALCULATED BASED ON THE SS INFILTRATION RATE ',
     +                 'IF STARTING WITH SS.'
         WRITE(iout,*)
      END IF
      lloc = llocsave
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'NOSURFLEAK') THEN
         Iseepsupress = 1
         WRITE(iout,*)
         WRITE(IOUT,'(A)')' SURFACE LEAKAGE WILL NOT BE SIMULATED '
         WRITE(iout,*)
         llocsave = lloc
      END IF
! lloc = llocsave
! CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
! IF(LINE(ISTART:ISTOP).EQ.'ONLYET') THEN
! IETBUD = 1
! WRITE(iout,*)
! WRITE(IOUT,'(A)')' ONLY ET WILL BE WRITTEN TO THE UNFORMATTED 
!+                      BUDGET FILES '
! WRITE(iout,*)
! END IF
      IF ( ITHTIFLG.GT.0 .OR. ITHTRFLG.GT.0 .OR. Iseepsupress.GT.0 ) 
     +     CALL URDCOM(In, IOUT, line)
      lloc = 1
      CALL URWORD(line, lloc, istart, istop, 2, NUZTOP, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, IUZFOPT, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, IRUNFLG, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, IETFLG, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, IUZFCB1, r, IOUT, In)
      CALL URWORD(line, lloc, istart, istop, 2, IUZFCB2, r, IOUT, In)
      IUZFB22 = IUZFCB2
      IUZFB11 = IUZFCB1
      NTRAIL = 1
      NSETS = 1
C
C2------READ UNSATURATED FLOW FLAGS WHEN IUZFOPT IS GREATER THAN ZERO.
      IF ( IUZFOPT.GT.0 ) THEN
        CALL URWORD(line, lloc, istart, istop, 2, NTRAIL, r, IOUT, In)
        CALL URWORD(line, lloc, istart, istop, 2, NSETS, r, IOUT, In)
      END IF
      CALL URWORD(line, lloc, istart, istop, 2, NUZGAG, r, IOUT, In)
      i=1
      CALL URWORD(line, lloc, istart, istop, 3, i, SURFDEP, IOUT, In)
!      CALL URWORD(line, lloc, istart, istop, 2, Iseepsupress, r, IOUT, 
!     +            In)
      RTSOLUTE = 0
!      CALL URWORD(line, lloc, istart, istop, 2, RTSOLUTE, r, IOUT, In)
C
C3------CHECK FOR ERRORS.
      IF ( SURFDEP.LT.ZEROD6 ) SURFDEP = ZEROD6
      IF ( IUZFOPT.GT.2 ) THEN
        WRITE (IOUT, 9002)
 9002   FORMAT (//' VERTICAL FLOW IN VADOSE ZONE IS ', 
     +          'ACTIVE AND IUZFOPT IS NOT 1 OR 2 '
     +          //, ' PLEASE CHECK INPUT  --  SETTING UZF PACKAGE TO ', 
     +          'INACTIVE'///)
        In = 0
        RETURN
      ELSE IF ( IUZFOPT.LE.0 ) THEN
        WRITE (IOUT, 9003)
 9003   FORMAT (//' UNSATURATED FLOW IN VADOSE ZONE IS IGNORED ', 
     +          //'RECHARGE TO GROUND WATER IS EQUAL TO SPECIFIED ', 
     +          'INFILTRATION RATE MINUS REJECTED RECHARGE'///)
      WRITE (IOUT, 9999)
 9999   FORMAT (//'***WARNING*** IUZFOPT IS ZERO. UNSATURATED  ',/ 
     +          'VARIABLES, EPS, THTS, and THTI ARE NOT READ. ',/ 
     +          'IF ET IS ACTIVE EXTWC IS NOT READ.'///)
      END IF
      IF ( ABS(IUZFOPT).EQ.2 .AND. Iunitbcf.GT.0 ) THEN
        WRITE (IOUT, 9004) IUZFOPT
 9004   FORMAT (//'BCF PACKAGE IS ACTIVE AND IUZFOPT = ', I5, 
     +          ' -- '//' ABSOLUTE VALUE OF IUZFOPT MUST EQUAL 1 ', 
     +          'FOR BCF  --  SETTING UZF PACKAGE TO INACTIVE'///)
        In = 0
        RETURN
      END IF
      IF ( NTRAIL.LT.0 ) THEN
        WRITE (IOUT, 9005)
 9005   FORMAT (//' NUMBER OF TRAILING WAVES IS LESS THAN ZERO'//
     +          ' --SETTING NTRAIL TO A POSITIVE VALUE'///)
        NTRAIL = ABS(NTRAIL)
      END IF
      IF ( NTRAIL.EQ.0 ) THEN
        WRITE (IOUT, 9006)
 9006   FORMAT (//' VERTICAL FLOW IN VADOSE ZONE IS ', 
     +          'ACTIVE AND NUMBER OF TRAILING WAVES IS ZERO  -- '//
     +          ' PLEASE CHECK INPUT  --  SETTING UZF PACKAGE TO ', 
     +          'INACTIVE'///)
        In = 0
        RETURN
      END IF
      IF ( IUZFOPT.GT.0 .AND. NSETS.LT.20 ) THEN
        WRITE (IOUT, 9007)
 9007   FORMAT (//' VERTICAL FLOW THROUGH UNSATURATED ZONE IS ', 
     +          'ACTIVE AND NUMBER OF WAVE SETS IS LESS THAN 20-- ', 
     +          ' RESETTING THE NUMBER OF WAVE SETS TO BE 20'///)
        NSETS = 20
      END IF
!      IF ( ABS(IUZFOPT).EQ.2 .AND. Iunithuf.GT.0 ) THEN
!        WRITE (IOUT, 9008) IUZFOPT
! 9008   FORMAT (//' VERTICAL FLOW IN VADOSE ZONE IS ', 
!     +          'ACTIVE, HUF PACKAGE IS ACTIVE AND IUZFOPT = ', I5, 
!     +          ' -- '//' ABSOLUTE VALUE OF IUZFOPT MUST EQUAL 1 ', 
!     +          'FOR HUF  --  SETTING UZF PACKAGE TO INACTIVE'///)
!        In = 0
!        RETURN
!      END IF
C
C4------ALLOCATE SPACE FOR UNSATURATED FLOW.
      ALLOCATE (IUZFBND(NCOL,NROW))
      IUZFBND = 0
C
C10-----READ IN BOUNDARY ARRAY FOR UNSATURATED FLOW. 
      CALL ReadArray(IUZFBND, aname(1), NROW, NCOL, 0, In, IOUT)
!     ALLOCATE ONLY CELLS THAT HAVE A NON-ZERO VALUE FOR IUZFBND
      NUMACTIVE = 0
      DO ir = 1, NROW
        DO ic = 1, NCOL
          IF ( IUZFBND(ic, ir).GT.0 ) NUMACTIVE = NUMACTIVE + 1
        END DO
      END DO
      IF ( IUZFOPT.EQ.1 .OR. IUZFOPT.EQ.2 ) THEN
        IUZM = NUMACTIVE
        NWAV = NTRAIL*(NSETS+1)
        NUZCL = NCOL
        NUZRW = NROW
C
C5------SET DIMENSIONS FOR UNSATURATED ZONE ARRAYS TO 1 IF NO
C      UNSATURATED ZONE
      ELSE
        IUZM = 1
        NWAV = 1
        NUZCL = 1
        NUZRW = 1
      END IF
      ALLOCATE (CHECKTIME(NWAV), MORE(NWAV))
C6------CALCULATE SPACE USED FOR LISTING UNSATURATED MOISTURE PROFILES.
      IF ( NUZGAG.GT.0 ) THEN
        NUZGAGAR = NUZGAG
      ELSE
        NUZGAGAR = 1
      END IF
C
C7------ALLOCATE SPACE FOR ARRAYS AND INITIALIZE.
      ALLOCATE (VKS(NCOL,NROW))
      VKS = 0.0 
      ALLOCATE (EPS(NUZCL,NUZRW), THTS(NUZCL,NUZRW), THTI(NUZCL,NUZRW))
      EPS = 0.0
      THTS = 0.0
      THTI = 0.0
      TOTRUNOFF = 0.0
      ALLOCATE (THTR(NUZCL,NUZRW))
      THTR = 0.0
      ALLOCATE (FINF(NCOL,NROW),PETRATE(NCOL,NROW),UZFETOUT(NCOL,NROW))
      ALLOCATE (GWET(NCOL,NROW), FNETEXFIL(NCOL,NROW))
      IF ( IETBUD.GT.0 ) THEN
        ALLOCATE (CUMGWET(NCOL,NROW))
      ELSE
        ALLOCATE (CUMGWET(1,1))
      END IF
      FINF = 0.0
      PETRATE = 0.0
      UZFETOUT = 0.0
      GWET = 0.0
      CUMGWET = 0.0
      FNETEXFIL = 0.0
      ALLOCATE (FBINS(52))
      FBINS = 0.0
      ALLOCATE (ROOTDPTH(NCOL,NROW))
      ROOTDPTH = 0.0
      ALLOCATE (WCWILT(NUZCL,NUZRW))
      WCWILT = 0.0
      ALLOCATE (SEEPOUT(NCOL,NROW), EXCESPP(NCOL,NROW))
      IF ( RTSOLUTE.GT.0 ) THEN
        ALLOCATE (AIR_ENTRY(NCOL,NROW), H_ROOT(NCOL,NROW))
      ELSE
        ALLOCATE (AIR_ENTRY(1,1), H_ROOT(1,1))
      END IF
      ALLOCATE (REJ_INF(NCOL,NROW))
      SEEPOUT = 0.0
      EXCESPP = 0.0
      REJ_INF = 0.0
      AIR_ENTRY = -16.0
      H_ROOT = -15000.0
      ALLOCATE (IUZLIST(4, NUZGAGAR))
      IUZLIST = 0
      ALLOCATE (NWAVST(NUZCL,NUZRW))
      NWAVST = 1
      ALLOCATE (CUMUZVOL(5))
      CUMUZVOL = 0.0D0
      ALLOCATE (UZTSRAT(7))
      UZTSRAT = 0.0D0
      ALLOCATE (UZTOTBAL(NCOL,NROW,7))
      UZTOTBAL = 0.0D0
crgn changed allocation 10/23/06
      ALLOCATE (UZFLWT(NCOL,NROW))
      UZFLWT = 0.0D0
      ALLOCATE (UZSTOR(NUZCL,NUZRW))
      UZSTOR = 0.0D0
      ALLOCATE (DELSTOR(NUZCL,NUZRW))
      DELSTOR = 0.0D0
cdep changed allocation of UZOLSFLX 7/30/08
      ALLOCATE (UZOLSFLX(NCOL,NROW))
      UZOLSFLX = 0.0D0
      ALLOCATE (HLDUZF(NCOL,NROW))
      HLDUZF = 0.0D0
      ALLOCATE (IUZHOLD(2, NCOL*NROW))
      nrnc = 1
      DO irhld = 1, NROW
        DO ichld = 1, NCOL
          IUZHOLD(1, nrnc) = irhld
          IUZHOLD(2, nrnc) = ichld
          nrnc = nrnc + 1
        END DO
      END DO
      ALLOCATE (ITRLSTH(NWAV))
      ITRLSTH = 0
      ALLOCATE (UZDPIT(NWAV))
      UZDPIT = 0.0D0
      ALLOCATE (UZDPST(NWAV,IUZM))
      UZDPST = 0.0D0
      ALLOCATE (UZTHIT(NWAV))
      UZTHIT = 0.0D0
      ALLOCATE (UZTHST(NWAV,IUZM))
      UZTHST = 0.0D0
      ALLOCATE (UZSPIT(NWAV))
      UZSPIT = 0.0D0
      ALLOCATE (UZSPST(NWAV,IUZM))
      UZSPST = 0.0D0
      ALLOCATE (UZFLIT(NWAV))
      UZFLIT = 0.0D0
      ALLOCATE (UZFLST(NWAV,IUZM))
      UZFLST = 0.0
      ALLOCATE (LTRLIT(NWAV))
      LTRLIT = 0
      ALLOCATE (LTRLST(NWAV,IUZM))
      LTRLST = 0
      ALLOCATE (ITRLIT(NWAV))
      ITRLIT = 0
      ALLOCATE (ITRLST(NWAV,IUZM))
      ITRLST = 0
      IF ( RTSOLUTE.GT.0 ) THEN
        ALLOCATE (RTSOLWC(NLAY,NCOL*NROW))
        RTSOLWC = 0
        ALLOCATE (RTSOLFL(NLAY,NCOL*NROW))
        RTSOLFL = 0
        ALLOCATE (RTSOLDS(NLAY,NCOL*NROW))
        RTSOLDS = 0
      ELSE
        ALLOCATE (RTSOLWC(1,1))
        RTSOLWC = 0
        ALLOCATE (RTSOLFL(1,1))
        RTSOLFL = 0
        ALLOCATE (RTSOLDS(1,1))
        RTSOLDS = 0
      END IF
      
C
C8------PRINT OPTION CODE WHEN NUZTOP IS WITHIN SPECIFIED RANGE.
      IF ( IUZFOPT.LE.0 ) THEN
        iuzflg = 0
      ELSE
        iuzflg = 1
      END IF
C
C8b-----Set flag for determining if FINF will be provided by PRMS.
C       A value of zero means that FINF will not be set by PRMS.
      IGSFLOW = 0
      
      IF ( NUZTOP.GE.1 .AND. NUZTOP.LE.3 ) THEN
        IF ( NUZTOP.EQ.1 ) WRITE (IOUT, 9009)
 9009   FORMAT (' OPTION 1 -- RECHARGE IN UZF TO TOP LAYER ONLY ')
        IF ( NUZTOP.EQ.2 ) WRITE (IOUT, 9010)
 9010   FORMAT (' OPTION 2 -- RECHARGE IN UZF TO SPECIFIED NODE ', 
     +          'IN EACH VERTICAL COLUMN')
        IF ( NUZTOP.EQ.3 ) WRITE (IOUT, 9011)
 9011   FORMAT (' OPTION 3 -- RECHARGE IN UZF TO HIGHEST ACTIVE ', 
     +          'NODE IN EACH VERTICAL COLUMN')
C
C9------STOP SIMULATION IF NUZTOP IS NOT WITHIN SPECIFIED RANGE.
      ELSE
        WRITE (IOUT, 9012) NUZTOP
 9012   FORMAT (1X, 'ILLEGAL RECHARGE OPTION CODE IN UZF (NUZTOP = ', 
     +          I5, ') -- SIMULATION ABORTING')
        CALL USTOP(' ')
      END IF
C
C11-----READ STREAM AND LAKE ARRAY FOR ROUTING OVERLAND FLOW.
      IF ( IRUNFLG.GT.0 ) THEN
        ALLOCATE(IRUNBND(NCOL,NROW)) 
        IRUNBND = 0
        CALL ReadArray(IRUNBND, aname(2), NROW, NCOL, 0, 
     +                                In, IOUT)
      ELSE
        ALLOCATE(IRUNBND(1,1)) 
        IRUNBND = 0
      END IF
C
C5B------READ AND SET VALUES FOR SOLUTE ROUTING IN UNSATURATED ZONE
C
      IF ( RTSOLUTE.GT.0 ) THEN
        ALLOCATE(GRIDSTOR(NCOL,NROW,NLAY))   
        ALLOCATE(GRIDET(NCOL,NROW,NLAY))
        GRIDSTOR = 0.0D0  
        GRIDET = 0.0D0
      ELSE
        ALLOCATE(GRIDSTOR(1,1,1),GRIDET(1,1,1))
        GRIDSTOR(1,1,1) = 0.0
        GRIDET(1,1,1) = 0
      END IF
C
C12-----READ VERTICAL HYDRAULIC CONDUCTIVITY FROM UZF INPUT FILE.
      IF ( IUZFOPT.EQ.1 .OR. IUZFOPT.LE.0 ) THEN
        CALL ReadArray(VKS, aname(6), NROW, NCOL, 0, In, IOUT)
C
C13-----CHECK FOR ERRORS IN VERTICAL HYDRAULIC CONDUCTIVITY
        DO nrck = 1, NROW
          DO ncck = 1, NCOL
            iflgbnd = 1
            IF ( IUZFBND(ncck, nrck).NE.0 ) THEN
              IF ( VKS(ncck, nrck).LT.CLOSEZERO ) THEN
                WRITE (IOUT, 9013) nrck, ncck
 9013           FORMAT (1X/, 'SATURATED VERTICAL K FOR CELL AT ROW ', 
     +                  I5, ', COL. ', I5, ' IS LESS THAN OR EQUAL TO ',
     +                  'ZERO-- SETTING UNSATURATED FLOW IN CELL TO ', 
     +                  'INACTIVE')
                iflgbnd = 0
              END IF
            END IF
            IF ( iflgbnd.EQ.0 ) IUZFBND(ncck, nrck) = 0
          END DO
        END DO
      END IF
      IF ( iuzflg.GT.0 ) THEN
C
C14-----READ BROOKS-COREY EPSILON ASSUMING IT IS CONSTANT THROUGHOUT
C        VERTICAL COLUMN.
        CALL ReadArray(EPS, aname(5), NUZRW, NUZCL, 0, In, IOUT)
C
C15-----READ SATURATED WATER CONTENT FOR UNSATURATED ZONE ASSUMING IT
C       IS CONSTANT THROUGHOUT VERTICAL COLUMN.
        CALL ReadArray(THTS, aname(3), NUZRW, NUZCL, 0, In, IOUT)
C
        IF ( ITHTRFLG.GT.0 )THEN
          CALL ReadArray(THTR, aname(8), NUZRW, NUZCL, 0, In, IOUT)
        END IF
C
C15-----READ AIR ENTRY PRESSURE FOR UNSATURATED ZONE ASSUMING IT
C       IS CONSTANT THROUGHOUT VERTICAL COLUMN.
!        CALL U2DREL(AIR_ENTRY, aname(3), NUZRW, NUZCL, 0, In, IOUT)
C
C
C15-----READ AIR ENTRY PRESSURE FOR UNSATURATED ZONE ASSUMING IT
C       IS CONSTANT THROUGHOUT VERTICAL COLUMN.
!        CALL U2DREL(H_ROOT, aname(3), NUZRW, NUZCL, 0, In, IOUT)
C
C16-----READ INITIAL WATER CONTENT FOR UNSATURATED ZONE ASSUMING IT
C         IS CONSTANT THROUGHOUT VERTICAL COLUMN. DO NOT READ
C         INITIAL WATER CONTENT IF PERIOD IS STEADY STATE.
        IF ( ISSFLG(1).EQ.0 .OR. ITHTIFLG.GT.0 )  
     +       CALL ReadArray(THTI, aname(4), NUZRW, NUZCL, 0, In, IOUT)
C
C17-----CHECK FOR ERRORS IN EPS, THTS, AND THTI ARRAYS.
        DO nrck = 1, NUZRW
          DO ncck = 1, NUZCL
            iflgbnd = 1
            IF ( IUZFBND(ncck, nrck).GT.0 ) THEN            
              IF ( THTS(ncck, nrck).LT.CLOSEZERO ) THEN
                WRITE (IOUT, 9014) nrck, ncck, THTS(ncck, nrck)
 9014           FORMAT (1X/, 'SATURATED WATER CONTENT FOR CELL ', 
     +                  'AT ROW ', I5, ', COL. ', I5, 
     +                  ' IS LESS THAN OR EQUAL', 
     +                  ' TO ZERO-- SETTING UNSATURATED FLOW IN ', 
     +                  'CELL TO INACTIVE', F12.3)
                iflgbnd = 0
              END IF
              IF ( ISSFLG(1).EQ.0 ) THEN
                IF ( THTI(ncck, nrck).LT.CLOSEZERO ) THEN
                  WRITE (IOUT, 9015) nrck, ncck, THTI(ncck, nrck)
 9015             FORMAT (1X/, 'INITIAL WATER CONTENT FOR CELL AT ', 
     +                    'ROW ', I5, ', COL. ', I5, 
     +                    ' IS LESS THAN OR EQUAL', 
     +                    ' TO ZERO-- SETTING UNSATURATED FLOW IN CELL '
     +                    , 'TO INACTIVE', F12.3)
                  iflgbnd = 0
                END IF
              END IF
              IF ( EPS(ncck, nrck).LT.CLOSEZERO ) THEN
                WRITE (IOUT, 9016) nrck, ncck, EPS(ncck, nrck)
 9016           FORMAT (1X/, 'BROOKS-COREY EPSILON FOR CELL AT ROW ', 
     +                  I5, ', COL. ', I5, ' IS LESS THAN OR EQUAL TO ',
     +                  'ZERO-- SETTING UNSATURATED FLOW IN CELL TO ', 
     +                  'INACTIVE', F12.3)
                iflgbnd = 0
              END IF
!              IF ( AIR_ENTRY(ncck, nrck).GT.-CLOSEZERO ) THEN
!                WRITE (IOUT, 9029) nrck, ncck, AIR_ENTRY(ncck, nrck)
! 9029           FORMAT (1X/, 'AIR ENTERY PRESSURE FOR CELL AT ROW ', 
!     +              I5, ', COL. ', I5, ' IS GREATER THAN OR EQUAL TO ',
!     +              'ZERO-- SETTING UNSATURATED FLOW IN CELL TO ', 
!     +                  'INACTIVE', F12.3)
!                iflgbnd = 0
!              END IF
!              IF ( H_ROOT(ncck, nrck).GT.-CLOSEZERO ) THEN
!                WRITE (IOUT, 9030) nrck, ncck, H_ROOT(ncck, nrck)
! 9030           FORMAT (1X/, 'AIR ENTERY PRESSURE FOR CELL AT ROW ', 
!     +              I5, ', COL. ', I5, ' IS GREATER THAN OR EQUAL TO ',
!     +              'ZERO-- SETTING UNSATURATED FLOW IN CELL TO ', 
!     +                  'INACTIVE', F12.3)
!                iflgbnd = 0
!              END IF
            END IF
            IF ( iflgbnd.EQ.0 ) IUZFBND(ncck, nrck) = 0
          END DO
        END DO
C
C18-----COMPUTE RESIDUAL WATER CONTENT (THTR) FOR ALL UNSATURATED CELLS
C       AS THE DIFFERENCE BETWEEN SATURATED WATER CONTENT AND SPECIFIC
C       YIELD OF UPPERMOST ACTIVE LAYER.
        IF ( IUZFOPT.GT.0 .AND. ITRSS.NE.0 .AND. ITHTRFLG.EQ.0 ) THEN
          DO nrth = 1, NROW
            isyflg = 1
            DO ncth = 1, NCOL 
              iuzflg = IUZFBND(ncth, nrth)
              IF ( iuzflg.GT.0 ) THEN
                nlth = 1
                ibndflg = 0
                DO WHILE ( ibndflg.EQ.0 )
                  ibndflg = IBOUND(ncth, nrth, nlth)
                  IF ( ibndflg.LT.1 ) nlth = nlth + 1
                  IF ( nlth.GT.NLAY ) ibndflg = -1
                END DO
C
C19-----SPECIFIC YIELD IS STORAGE CAPACITY DIVIDED BY AREA OF MODEL CELL.
                IF ( ibndflg.GT.0 ) THEN                        
                  sy = 0.0
                  IF ( Iunitlpf.GT.0 ) THEN
                    IF ( LAYHDT(nlth).GT.0 ) THEN
C use LPF SC2, Iunitlpf>0
                      sy = SCLPF(ncth, nrth, nlth)/
     +                     (DELR(ncth)*DELC(nrth))
                    ELSE
                      WRITE (IOUT, 9017) nlth, ncth, nrth 
9017                  FORMAT(1X,'PROGRAM TERMINATED-LAYTYP IN LPF '
     +                       , 'PACKAGE MUST BE GREATER THAN ZERO FOR '
     +                       , 'UPPERMOST ACTIVE CELL',/1X
     +                       , 'SO SPECIFIC YIELD CAN BE USED FOR '
     +                       , 'COMPUTING RESIDUAL WATER CONTENT-- '
     +                       , 'CELL LAYER,ROW,COLUMN: ',3I5) 
                     CALL USTOP(' ')
                    END IF
                  ELSE IF ( Iunitupw.GT.0 ) THEN
                    IF ( LAYHDT(nlth).EQ.0 ) THEN
                      THICK = BOTM(ncth, nrth, LBOTM(nlth)-1)-
     +                      BOTM(ncth, nrth, LBOTM(nlth))
                      sy = SC2UPW(ncth, nrth, nlth)/
     +                     (THICK*DELR(ncth)*DELC(nrth)) 
                    ELSE
                      sy = SC2UPW(ncth, nrth, nlth)/
     +                     (DELR(ncth)*DELC(nrth)) 
                    END IF 
                  ELSE IF ( Iunitbcf.GT.0 ) THEN                
                    IF ( LAYCON(nlth).LT.2 ) THEN
                      sy = SC1(ncth, nrth, nlth)/
     +                     (DELR(ncth)*DELC(nrth))
                    ELSE
                      kkrch = 0
                      DO k = 1, nlth
                        IF(LAYCON(K).EQ.3 .OR. LAYCON(K).EQ.2)
     +                   kkrch = kkrch + 1
                      END DO
                      sy = SC2(ncth, nrth, kkrch)/
     +                     (DELR(ncth)*DELC(nrth))
                    END IF
*                  ELSE IF ( Iunithuf.GT.0 ) THEN 
*                    sy = SC2HUF(ncth, nrth)
                  END IF
                  IF ( sy.GT.0.0 ) THEN
                    isyflg = 1
                    THTR(ncth, nrth) = THTS(ncth, nrth) - sy
                  ELSE
                    isyflg = 0
                    THTR(ncth, nrth) = 0.0D0
                  END IF
                ELSE IF ( ibndflg.EQ.0 ) THEN
                  isyflg = 0
                END IF
                test = (THTI(ncth, nrth)-THTR(ncth, nrth))
                IF ( test.LT.CLOSEZERO ) THTI(ncth,nrth) = 
     +               THTR(ncth,nrth)
                test = (THTS(ncth, nrth)-THTR(ncth, nrth))
                IF ( test.LT.CLOSEZERO ) THEN
                  WRITE (IOUT, 9028) nrth, ncth
                  IUZFBND(ncth, nrth) = 0
                END IF
              ELSE
                THTR(ncth, nrth) = 0.0D0
              END IF
9028  FORMAT (1X/, 'SATURATED WATER CONTENT FOR UPPERMOST ACTIVE ', 
     +                  'CELL AT ROW ', I5, ', COL. ', I5, ' IS LESS ', 
     +                  'THAN OR EQUAL TO RESIDUAL WATER CONTENT-- ',
     +                  'SETTING UNSATURATED FLOW IN CELL TO INACTIVE')
C
C20-----IF SPECIFIC YIELD IS 0 FOR UPPERMOST ACTIVE CELL AT COLUMN J
C       AND ROW I OR IF ALL LAYERS AT CELL ARE INACTIVE, SET
C       UNSATURATED FLOW AT CELL INACTIVE.
              IF ( isyflg.EQ.0 ) THEN
                WRITE (IOUT, 9018) nrth, ncth
 9018           FORMAT (1X/, 'SPECIFIC YIELD FOR UPPERMOST ACTIVE ', 
     +                  'CELL AT ROW ', I5, ', COL. ', I5, ' IS LESS ', 
     +                  'THAN OR EQUAL TO ZERO-- SETTING UNSATURATED ',
     +                  'FLOW IN CELL TO INACTIVE')
                IUZFBND(ncth, nrth) = 0
              END IF
            END DO
          END DO
        END IF   
      END IF
C
C21-----READ FILES FOR PRINTING TIME SERIES WATER CONTENT PROFILES.
      IF ( NUZGAG.GT.0 ) THEN
        WRITE (IOUT, 9019)
 9019   FORMAT (1X/, 'UNSATURATED FLOW RESULTS WILL BE PRINTED TO ', 
     +          'SEPARATE FILES FOR SELECTED MODEL CELLS, AND ', 
     +          'TIME STEPS AS DEFINED BY OUTPUT CONTROL PACKAGE ', 
     +          //, 'SELECTED MODEL CELLS ARE: '/' ROW NUMBER ', 
     +          ' COLUMN NUMBER  FORTRAN UNIT NUMBER  ', 
     +          'OUTPUT OPTION'//)
        igage = 1
        DO WHILE ( igage.LE.NUZGAG )
          READ (In, *) IUZLIST(1, igage) 
          IF( IUZLIST(1, igage) .GE. 0 ) THEN 
            BACKSPACE In
            READ (In, *) (IUZLIST(l, igage), l=1, 4)
          ELSE
            IUZLIST(3, igage) = -1*IUZLIST(1, igage)
            IUZLIST(4, igage) = 4  
            IUZLIST(1, igage) = 0
            IUZLIST(2, igage) = 0            
          END IF
          iuzrow = IUZLIST(1, igage)
          iuzcol = IUZLIST(2, igage)
          IF (iuzcol.GT.0 .AND. iuzrow.GT.0 ) THEN
            icheck = IUZFBND(iuzcol, iuzrow)
          ELSE
            icheck = 0
          END IF
          IF ( IUZLIST(4, igage).EQ.3 .AND. icheck.LE.0 ) THEN
            WRITE (IOUT,*) '**WARNING** Printing of water content ',
     +                     'profiles is not possible when IUZFOPT <= 0.'
            WRITE (IOUT,*) 'Removing gage',igage
            igage = igage - 1
            NUZGAG = NUZGAG - 1
          ELSE
            WRITE (IOUT, 9020) (IUZLIST(l, igage), l=1, 4)
 9020     FORMAT (1X, I7, 7X, I7, 12X, I7, 11X, I7)
C
C22-----DETERMINE IF ROW AND COLUMN NUMBERS ARE IN ACTIVE AREA OF
C        UNSATURATED FLOW.
            IF ( IUZLIST(4, igage) .LT. 4 ) THEN
              IF ( iuzrow.LT.1 .OR. iuzrow.GT.NROW) THEN
                WRITE (IOUT, 9021) iuzrow, igage
 9021         FORMAT (1X/, 'WARNING--- ROW NUMBER ', I7,  
     +              ' FOR RECORD ',I7, ' IS NOT A VALID NUMBER', 
     +              ' NO OUTPUT WILL BE PRINTED TO THIS FILE'/)
                IUZLIST(3, igage) = 0
              END IF
              IF ( iuzcol.LT.1 .OR. iuzcol.GT.NCOL ) THEN
                WRITE (IOUT, 9022) iuzcol, igage
 9022         FORMAT (1X/, 'WARNING--- COLUMN NUMBER ', I7, 
     +              ' FOR RECORD ', I7, ' IS NOT A VALID NUMBER', 
     +              ' NO OUTPUT WILL BE PRINTED TO THIS FILE'/)
                IUZLIST(3, igage) = 0
              END IF
              IF ( icheck.EQ.0 ) THEN
                WRITE (IOUT, 9023) igage
 9023         FORMAT (1X/, 'WARNING--- RECORD ', I7, ' IS NOT IN ', 
     +              'ACTIVE AREA OF UNSATURATED FLOW; ', 
     +              ' NO OUTPUT WILL BE PRINTED TO THIS FILE'/)
              IUZLIST(3, igage) = 0
              END IF
            END IF
          END IF
          igage = igage + 1
        END DO
C
C23-----WRITE HEADER FILES FOR CELLS WITH SELECTED OUTPUT.
        DO igage = 1, NUZGAG
          iuzrow = IUZLIST(1, igage)
          iuzcol = IUZLIST(2, igage)
          igunit = IUZLIST(3, igage)
          iuzopt = IUZLIST(4, igage)
          IF ( iuzcol.GT.0 ) iuzlay = IUZFBND(iuzcol, iuzrow)

! Disable writing of gage files for now
*          IF ( igunit.GT.0 ) THEN
*C
*C24----GET VARIABLE OUTTYPE.
*            SELECT CASE (iuzopt)
*C
*C25-----PRINT HEADER WHEN WRITING ONLY VOLUMES.
*            CASE (1)
*              WRITE (igunit, 9024) igage, iuzrow, iuzcol, iuzlay
*C
*C26-----PRINT HEADER WHEN WRITING VOLUMES AND RATES.
*            CASE (2)
*              WRITE (igunit, 9025) igage, iuzrow, iuzcol, iuzlay
*C
*C27-----PRINT HEADER FOR UNSATURATED-ZONE MOISTURE PROFILES
*            CASE (3)
*              WRITE (igunit, 9026) igage, iuzrow, iuzcol, iuzlay
*            CASE (4)
*              WRITE (igunit, 9027)
*            END SELECT
*          END IF

        END DO
      END IF
C
C28-----FORMATS
 9024 FORMAT (1X,'"LOCATION OF SPECIFIED CELL FOR PRINTING VOLUMES ', 
     +        'IN UNSATURATED ZONE: GAGE ', I4, ' ROW, COLUMN ', I4, 
     +        ',', I4, ' INITIAL LAYER ASSIGMENT ', I4,'"', /1X,
     +        '"DATA:  LAYER             TIME        GW-HEAD   ', 
     +        'UZ-THICKNESS CUM.-APL.-INF.   CUM.-INFILT.     CUM.-', 
     +        'RECH.    TOTAL-STOR.   STOR.-CHANGE    SURF.-LEAK. "')
 9025 FORMAT (1X,'"LOCATION OF SPECIFIED CELL FOR PRINTING VOLUMES ', 
     +      'AND RATES IN UNSATURATED ZONE: GAGE ', I4, 
     +      ' ROW, COLUMN ', I4, ',', I4, ' INITIAL LAYER ASSIGMENT ',
     +      I4,'"', /1X, '"DATA:  LAYER             TIME        ',
     +      'GW-HEAD   UZ-THICKNESS CUM.-APL.-INF.   CUM.-INFILT.   ', 
     +      '  CUM.-RECH.    TOTAL-STOR.   STOR.-CHANGE    SURF.-LEAK.',
     +      ' APL.-INF.-RATE    INFIL.-RATE     RECH.-RATE     ',
     +      'STOR.-RATE     SEEP.-RATE "', /)
 9026 FORMAT (1X,'"LOCATION OF SPECIFIED CELL FOR PRINTING VOLUMES ', 
     +        'IN UNSATURATED ZONE: GAGE ', I4, ' ROW, COLUMN ', I4, 
     +        ',', I4, ' INITIAL LAYER ASSIGMENT ', I4,'"', /1X, 
     +        '"DATA:  LAYER             TIME        GW-HEAD   ',
     +        'UZ-THICKNESS          DEPTH    WATER-CONT. "', /)
 9027 FORMAT (1X,'"UNSATURATED MASS BALANCE COMPONENTS FOR ENTIRE ',  
     +        'MODEL "',/1X,'"DATA:            TIME   APPLIED-INFIL.',
     +        '          RUNOFF   ACTUAL-INFIL.   SURFACE-LEAK.',
     +        '           UZ-ET           GW-ET',
     +        '   UZSTOR-CHANGE        RECHARGE "',/)
C
C29-----PRINT WARINING WHEN UNITS ARE UNDEFINED IN MODFLOW.
      IF ( ITMUNI.EQ.0 .OR. LENUNI .EQ. 0 ) THEN
        WRITE(IOUT,*)'****Units are undefined. This may cause ',
     +  'unfortunate results when using GSFLOW****'
      END IF
C
      IF ( Iunitlpf.GT.0 .OR. Iunitupw.GT.0 ) THEN 
        IF ( ABS(IUZFOPT).NE.1 ) CALL SGWF2UZF1VKS(Iunitlpf,
     +                                             Iunitupw)
      END IF
      fkmax = 86400.0
      fkmin = 1.0E-4
      IF ( ITMUNI.EQ.1 ) THEN
        fkmax = fkmax/86400.0
        fkmin = fkmin/86400.0
      ELSE IF ( ITMUNI.EQ.2 ) THEN
        fkmax = fkmax/1440.0
        fkmin = fkmin/1440.0
      ELSE IF ( ITMUNI.EQ.3 ) THEN
        fkmax = fkmax/24.0  
        fkmin = fkmin/24.0
      ELSE IF ( ITMUNI.EQ.5 ) THEN
        fkmax = fkmax*365.0
        fkmin = fkmin*365.0
      END IF    
      IF ( LENUNI.EQ.1 ) THEN
        fkmax = fkmax/0.3048
        fkmin = fkmin/0.3048
      ELSE IF ( LENUNI.EQ.3 ) THEN
        fkmax = fkmax*100.0
        fkmin = fkmin*100.0
      END IF 
      range=LOG(fkmax)-LOG(fkmin)
      finc=range/50.0
      FBINS(1) = LOG(fkmin)
C
C30-----BIN INFILTRATION RATES.
      DO ivol = 2, 51
        FBINS(ivol) = FBINS(ivol-1)+ finc
        FBINS(ivol-1) = EXP(FBINS(ivol-1))
      END DO
      FBINS(51) = EXP(FBINS(51))
C
C31-----SAVE POINTERS FOR GRID AND RETURN.
      CALL SGWF2UZF1PSV(Igrid)
      RETURN
      END SUBROUTINE GWF2UZF1AR
C
C------SUBROUTINE SGWF2UZF1VKS
      SUBROUTINE SGWF2UZF1VKS(Iunitlpf, Iunitupw) 
C     ******************************************************************
C     ASSIGN SATURATED VERTICAL HYDRAULIC CONDUCTIVITY ARRAY 
C     (VKS) IN UZF TO EQUAL VERTICAL HYDRAULIC CONDUCTIVITY IN LAYER-
C     PROPERTY FLOW PACKAGE
!--------REVISED FOR MODFLOW-2005 RELEASE 1.9, FEBRUARY 6, 2012
C     VERSION 1.0.5:  April 5, 2012
C     ******************************************************************
      USE GWFUZFMODULE, ONLY: VKS, IUZFBND, CLOSEZERO, NUZTOP
      USE GLOBAL,       ONLY: NCOL, NROW, NLAY, IOUT, IBOUND, BOTM, 
     +                        LAYHDT
      USE GWFLPFMODULE, ONLY: LAYVKA, VKA, HK
*      USE GWFHUFMODULE, ONLY: HGUVANI, NHUF, HKHUF=>HK, VKAH
      USE GWFLPFMODULE, ONLY: SCLPF=>SC2
      USE GWFUPWMODULE, ONLY: LAYVKAUPW, VKAUPW, HKUPW
      IMPLICIT NONE
C    ------------------------------------------------------------------
C    SPECIFICATIONS:
C    ------------------------------------------------------------------
C     ARGUMENTS
C    ------------------------------------------------------------------
      INTEGER Iunitlpf, Iunitupw
C    ------------------------------------------------------------------
C    LOCAL VARIABLES
C    ------------------------------------------------------------------
      INTEGER krck, ncck, nrck, iflgbnd, il, ill, nlayp1
      REAL celthick
C    ******************************************************************
      nlayp1 = NLAY + 1
C
C1------CHECK TO SEE IF UPPERMOST ACTIVE CELL IS CONVERTABLE AND
C       SET VKS EQUAL TO VKALPF FOR CORRESPONDING MODEL CELL.
      krck = 1
      DO nrck = 1, NROW
        DO ncck = 1, NCOL
!   RGN 6/22/09. Add coded to find upper-most active layer
          il = 0
          IF ( NUZTOP.EQ.1 .OR. NUZTOP.EQ.2 ) THEN
            il = IUZFBND(ncck, nrck)
            IF ( il.GT.0 ) THEN
              IF ( IBOUND(ncck, nrck, il).LT.1 ) il = 0
            END IF
          ELSE IF ( NUZTOP.EQ.3 ) THEN
            ill = 1
            il = 0
            DO WHILE ( ill.LT.nlayp1 )
              IF ( IBOUND(ncck, nrck, ill).GT.0 ) THEN
                il = ill
                EXIT
              ELSE IF ( IBOUND(ncck, nrck, ill).LT.0 ) THEN
                EXIT
              END IF
CRGN made il = 0 when all layers for column are inactive 2/21/08
              ill = ill + 1
            END DO
          END IF
          krck = il
          IF ( krck.NE.0 ) THEN
            IF ( IBOUND(ncck, nrck, krck).GT.0 ) THEN
              IF ( Iunitlpf.GT.0 ) THEN
                IF ( LAYHDT(krck).LT.1 ) THEN
                  WRITE (IOUT, *) 
     +                       'PROGRAM TERMINATED-LAYTYP MUST BE GREATER'
     +                       , ' THAN ZERO WHEN IUZFOPT IS 2'
                  CALL USTOP(' ')
                END IF
                IF ( LAYVKA(krck).EQ.0 ) THEN
                  VKS(ncck, nrck) = VKA(ncck, nrck, krck)
                ELSE
                  VKS(ncck, nrck) = HK(ncck, nrck, krck)/
     +                              VKA(ncck, nrck, krck)
                END IF
*              ELSE IF ( Iunithuf.GT.0 ) THEN
*                IF ( krck.GT.0 ) THEN
*                  celthick = BOTM(ncck, nrck, krck-1)-
*     +                       BOTM(ncck, nrck, krck)
*                END IF
*                IF ( HGUVANI(NHUF).LT.CLOSEZERO ) THEN
*                  IF ( celthick.GT.1.0e-7 ) THEN
*                    VKS(ncck, nrck) = VKAH(ncck, nrck, krck)/(celthick)
*                  END IF
*                ELSE
*                  VKS(ncck, nrck) = HKHUF(ncck, nrck, krck)/
*     +                              HGUVANI(NHUF)
*                END IF
              ELSEIF ( Iunitupw.GT.0 ) THEN
                IF ( LAYHDT(krck).GT.0 ) THEN
                  IF ( LAYVKAUPW(krck).EQ.0 ) THEN
                    VKS(ncck, nrck) = VKAUPW(ncck, nrck, krck)
                  ELSE
                    VKS(ncck, nrck) = HKUPW(ncck, nrck, krck)/
     +                                VKAUPW(ncck, nrck, krck)
                  END IF
                ELSE
                  WRITE (IOUT, *) 
     +                       'PROGRAM TERMINATED-LAYTYP MUST BE GREATER'
     +                       , ' THAN ZERO WHEN IUZFOPT IS 2'
                  CALL USTOP(' ')
                END IF
              END IF
              iflgbnd = 0
              IF ( IUZFBND(ncck, nrck).NE.0 ) THEN
              iflgbnd = 1
              IF ( VKS(ncck, nrck).LT.CLOSEZERO ) THEN
                WRITE (IOUT, 9013) nrck, ncck
 9013           FORMAT (1X/, 'SATURATED VERTICAL K FOR CELL AT ROW ', 
     +                  I5, ', COL. ', I5, ' IS LESS THAN OR EQUAL TO ',
     +                  'ZERO-- SETTING UNSATURATED FLOW IN CELL TO ', 
     +                  'INACTIVE')
                iflgbnd = 0
              END IF
            END IF
            IF ( iflgbnd.EQ.0 ) IUZFBND(ncck, nrck) = 0
            END IF
          END IF
        END DO
      END DO
C
C2------RETURN.
      RETURN
      END SUBROUTINE SGWF2UZF1VKS
C
C-------SUBROUTINE GWF2UZF1RP
      SUBROUTINE GWF2UZF1RP(In, Kkper, Iunitsfr, Igrid)
C     ******************************************************************
C     READ AND CHECK VARIABLES EACH STRESS PERIOD 
!--------REVISED FOR MODFLOW-2005 RELEASE 1.9, FEBRUARY 6, 2012
C     VERSION 1.0.5:  April 5, 2012
C     ******************************************************************
      USE GWFUZFMODULE
      USE GLOBAL,       ONLY: NCOL, NROW, NLAY, IOUT, ISSFLG, IBOUND, 
     +                        HNEW, DELR, DELC, BOTM, LBOTM
      USE GWFSFRMODULE,ONLY:RECHSAVE
      IMPLICIT NONE
C     -----------------------------------------------------------------
C     SPECIFICATIONS:
C     -----------------------------------------------------------------
C     ARGUMENTS
C     -----------------------------------------------------------------
      INTEGER In, Kkper, Iunitsfr, Igrid
C     -----------------------------------------------------------------
C     LOCAL VARIABLES
C     -----------------------------------------------------------------
      DOUBLE PRECISION h
      DOUBLE PRECISION thtrcell
      DOUBLE PRECISION bottom, celtop, slen, width, etdpth, surfinf
      DOUBLE PRECISION thick, surfpotet, top
      INTEGER ic, iflginit, il, ilay, ill, ir, iss, jk, l, ncck,
     +        nrck, nuzf, ll, uzlay, land, nlayp1
      CHARACTER(LEN=24) aname(4)
      DATA aname(1)/' AREAL INFILTRATION RATE'/
      DATA aname(2)/'                 ET RATE'/
      DATA aname(3)/'     ET EXTINCTION DEPTH'/
      DATA aname(4)/'EXTINCTION WATER CONTENT'/
C     -----------------------------------------------------------------
      nlayp1 = NLAY + 1
C
C1------SET POINTERS FOR THE CURRENT GRID.
      CALL SGWF2UZF1PNT(Igrid)
C      
C2------READ INFILTRATION RATES FOR UZF CELLS AT THE BEGINNING OF EACH
C       STRESS PERIOD.
      iss = ISSFLG(Kkper)
      iflginit = 0
      READ (In, *) nuzf
      IF ( nuzf.LT.0 ) THEN
        WRITE (IOUT, *) 'USING INFILTRATION RATE FROM PREVIOUS STRESS '
        WRITE (IOUT, *) 'PERIOD.', 'CURRENT PERIOD IS: ', Kkper
      ELSE
C
C3------READ IN ARRAY FOR INFILTRATION RATE.
        CALL ReadArray(FINF, aname(1), NROW, NCOL, 0, In, IOUT)
C
C4------CHECK FOR NEGATIVE INFILTRATION RATES.
        DO nrck = 1, NROW
          DO ncck = 1, NCOL
            IF ( IUZFBND(ncck, nrck).NE.0 ) THEN
              surfinf = FINF(ncck, nrck)
              IF ( FINF(ncck, nrck).LT.0.0 ) THEN
                WRITE (IOUT, 9002) nrck, ncck
 9002           FORMAT (1X/, 'INFILTRATION RATE FOR CELL AT ROW ', I5, 
     +                ', COLUMN ', I5, ' IS LESS THAN ZERO ', 
     +                'ZERO-- SETTING RATE TO ZERO ')
                FINF(ncck, nrck) = 0.0
C
C5------SET INFILTRATION RATE TO SATURATED VERTICAL K WHEN RATE IS
C        GREATER THAN K AND ROUTE EXCESS WATER TO STREAM IF 
C        IRUNFLG IS NOT EQUAL TO ZERO.
              ELSE IF ( FINF(ncck, nrck).GT.VKS(ncck, nrck) ) THEN
                EXCESPP(ncck, nrck) =  (FINF(ncck, nrck) - 
     +                      VKS(ncck, nrck))*DELC(nrck)*DELR(ncck)
                FINF(ncck, nrck) = VKS(ncck, nrck)
              END IF
            END IF
            IF ( IUNITSFR.GT.0 ) RECHSAVE(ncck, nrck) = FINF(ncck, nrck)
          END DO
        END DO
      END IF
      IF ( IETFLG.NE.0 ) THEN
        READ (In, *) nuzf
        IF ( nuzf.LT.0 ) THEN
          WRITE (IOUT, 9004) Kkper
 9004     FORMAT (1X/, 'USING ET RATE FROM PREVIOUS STRESS ', 
     +            'PERIOD. CURRENT PERIOD IS: ', I7)
        ELSE
C
C6------READ IN ARRAY FOR ET RATE.
          CALL ReadArray(PETRATE, aname(2), NROW, NCOL, 0, In, IOUT)
C
C7-----CHECK FOR NEGATIVE ET RATES.
          DO nrck = 1, NROW
            DO ncck = 1, NCOL
C8----ZERO ET ARRAY THAT IS PASSED TO PRMS.
              UZFETOUT(ncck, nrck) = 0.0
              IF ( IUZFBND(ncck, nrck).EQ.0 ) THEN
                surfpotet = PETRATE(ncck, nrck)
                IF ( surfpotet.GT.0.0 ) PETRATE(ncck, nrck) = 0.0
              ELSE IF ( PETRATE(ncck, nrck).LT.0.0 ) THEN
                WRITE (IOUT, 9005) nrck, ncck
 9005           FORMAT (1X/, 'POTENTIAL ET RATE FOR CELL AT ROW ', I5, 
     +                  ', COLUMN ', I5, ' IS LESS THAN ZERO ', 
     +                  'ZERO-- SETTING RATE TO ZERO ')
                PETRATE(ncck, nrck) = 0.0
              END IF
            END DO
          END DO
        END IF
        READ (In, *) nuzf
        IF ( nuzf.LT.0 ) THEN
          WRITE (IOUT, 9006) Kkper
 9006     FORMAT (/1x, 'USING ET EXTINCTION DEPTH FROM PREVIOUS ',
     +            'STRESS PERIOD. CURRENT PERIOD IS: ', I7)
        ELSE
C
C9------READ IN ARRAY FOR ET EXTINCTION DEPTH.
          CALL ReadArray(ROOTDPTH, aname(3), NROW, NCOL, 0, In, IOUT)
C
C10-----CHECK FOR NEGATIVE ET EXTINCTION DEPTH.
          DO nrck = 1, NROW
            DO ncck = 1, NCOL
              IF ( IUZFBND(ncck, nrck).EQ.0 ) THEN
                etdpth = ROOTDPTH(ncck, nrck)
                IF ( etdpth.GT.0.0 ) ROOTDPTH(ncck, nrck) = 0.0
              ELSE IF ( ROOTDPTH(ncck, nrck).LT.CLOSEZERO ) THEN
                WRITE (IOUT, 9007) nrck, ncck
 9007           FORMAT (1X/, 'ROOTING DEPTH FOR CELL AT ROW ', I5, 
     +                  ', COLUMN ', I5, ' IS LESS THAN OR EQUAL TO ', 
     +                  'ZERO-- SETTING DEPTH TO ZERO ')
!                ROOTDPTH(ncck, nrck) = 1.0
                ROOTDPTH(ncck, nrck) = 0.0
              END IF
!   RGN 6/22/09. Add coded to find upper-most active layer
              il = 0
              IF ( NUZTOP.EQ.1 .OR. NUZTOP.EQ.2 ) THEN
                il = IUZFBND(ncck, nrck)
                IF ( il.GT.0 ) THEN
                  IF ( IBOUND(ncck, nrck, il).LT.1 ) il = 0
                END IF
              ELSE IF ( NUZTOP.EQ.3 ) THEN
                ill = 1
                il = 0
                DO WHILE ( ill.LT.nlayp1 )
                  IF ( IBOUND(ncck, nrck, ill).GT.0 ) THEN
                    il = ill
                    EXIT
                  ELSE IF ( IBOUND(ncck, nrck, ill).LT.0 ) THEN
                    EXIT
                  END IF
CRGN made il = 0 when all layers for column are inactive 2/21/08
                  ill = ill + 1
                END DO
              END IF
              land = ABS(IUZFBND(ncck, nrck))
!
              IF ( il.GT.0 .AND. land.GT.0 ) THEN
                thick = BOTM(ncck, nrck,LBOTM(land)-1)-
     +                  BOTM(ncck, nrck,LBOTM(il))
                IF ( ROOTDPTH(ncck, nrck).GT.0.9*thick ) THEN
                  ROOTDPTH(ncck, nrck) = 0.9*thick
                  WRITE (IOUT, 222) nrck, ncck
                END IF
              END IF
            END DO
          END DO
        END IF
  222 FORMAT('***WARNING*** ET EXTINCTION DEPTH IS BELOW LAYER BOTTOM',/
     +       'RESETTING DEPTH TO 90% OF CELL THICKNESS FOR ROW ',I7,
     +       ' COLUMN ',I7)
C
C11-----SKIP READING OF EXTINCTION WATER CONTENT ARRAY WHEN 
C         IUZFOPT IS ZERO.
        IF ( IUZFOPT.GT.0 ) THEN
          READ (In, *) nuzf
          IF ( nuzf.LT.0 ) THEN
            WRITE (IOUT, 9012) Kkper
 9012       FORMAT (/1X, 'USING EXTINCTION WATER CONTENT FROM PREVIOUS',
     +              ' STRESS PERIOD. CURRENT PERIOD IS: ', I7)
          ELSE
C
C12-----READ IN ARRAY FOR ET EXTINCTION DEPTH.
            CALL ReadArray(WCWILT, aname(4), NROW, NCOL, 0, In, IOUT)
C
C13-----CHECK FOR EXTINCTION WATER CONTENT LESS THAN RESIDUAL WATER
C         CONTENT.
            DO nrck = 1, NUZRW
              DO ncck = 1, NUZCL
                IF ( IUZFBND(ncck, nrck).LT.1 .AND. WCWILT(ncck, nrck)
     +               .GT.0 ) WCWILT(ncck, nrck) = 0.0
                IF ( iss.NE.0 .AND. IUZFBND(ncck, nrck).GT.0 ) THEN
                  IF ( WCWILT(ncck, nrck).LT.THTR(ncck, nrck)+ZEROD9 )
     +                 THEN
                    WRITE (IOUT, 9014) nrck, ncck
 9014               FORMAT (1X/, 'EXTINCTION WATER CONTENT FOR ', 
     +                      'CELL AT ROW ', I5, ', COLUMN ', I5, 
     +                      ' IS LESS THAN RESIDUAL WATER CONTENT', 
     +                      '-- SETTING EXTINCTION WATER ', 
     +                      'CONTENT EQUAL TO RESIDUAL WATER ', 
     +                      'CONTENT')
                    WCWILT(ncck, nrck) = THTR(ncck, nrck) + ZEROD9
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF
C13B-----SEARCH FOR UPPERMOST ACTIVE CELL.
      DO ir = 1, NROW
        DO ic = 1, NCOL
          IF ( IUZFBND(ic, ir).NE.0 ) THEN
            il = 0
            IF ( NUZTOP.EQ.1 .OR. NUZTOP.EQ.2 ) THEN
              il = ABS(IUZFBND(ic, ir))
              IF ( il.GT.0 ) THEN
                IF ( IBOUND(ic, ir, il).LT.1 ) il = 0
              ELSE
                il = 0
              END IF
              IF ( IL.EQ.0 ) IUZFBND(ic, ir) = 0
            ELSE IF ( NUZTOP.EQ.3 ) THEN
              ill = 1
              il = 0
              DO WHILE ( ill.LT.nlayp1 )
                IF ( IBOUND(ic, ir, ill).GT.0 ) THEN
                  il = ill
                  EXIT
                ELSE IF ( IBOUND(ic, ir, ill).LT.0 ) THEN
                  EXIT
                END IF
                ill = ill + 1
              END DO
            END IF
          END IF
        END DO
      END DO
C
C14------INITIALIZE UNSATURATED ZONE IF ACTIVE.
C
C15------SET FLAGS FOR STEADY STATE OR TRANSIENT SIMULATIONS.
      IF ( Kkper.GT.2 ) THEN
        iflginit = 0
      ELSE IF ( Kkper.EQ.1 ) THEN
        iflginit = 1
      ELSE
        IF ( iss.EQ.0 .AND. ISSFLG(Kkper-1).NE.0 )
     +       iflginit = 2
      END IF
      IF ( iflginit.GE.1 ) THEN
        l = 0
        DO ll = 1, NUMCELLS
          ir = IUZHOLD(1, ll)
          ic = IUZHOLD(2, ll)
          IF ( IUZFBND(ic,ir).GT.0 ) THEN
            l = l + 1
C
C16-----SEARCH FOR UPPERMOST ACTIVE CELL.
            IF ( IUZFBND(ic, ir).GT.0 ) THEN
              il = 0 
              IF ( NUZTOP.EQ.1 .OR. NUZTOP.EQ.2 ) THEN
                il = IUZFBND(ic, ir)
                IF ( il.GT.0 ) THEN
                  IF ( IBOUND(ic, ir, il).LT.1 ) il = 0
                ELSE
                  il = 0
                END IF
              ELSE IF ( NUZTOP.EQ.3 ) THEN
                ill = 1
                il = 0
                DO WHILE ( ill.LT.nlayp1 )
                  IF ( IBOUND(ic, ir, ill).GT.0 ) THEN
                    il = ill
                    EXIT
                  ELSE IF ( IBOUND(ic, ir, ill).LT.0 ) THEN
                    EXIT
                  END IF
                  ill = ill + 1
                END DO
              END IF
C
C16B-----SEARCH FOR UPPER MOST ACTIVE CELL WITH A WATER LEVEL.
              ilay = il
              IF ( il.GT.0 ) THEN
                IF ( IBOUND(ic, ir, il).GT.0 ) THEN
                  TOPCELL: DO WHILE ( ilay.LT.nlayp1 )
!                    IF ( HNEW(ic, ir, ilay).LE.BOTM(ic,ir,ilay) ) THEN
!                      ilay = ilay + 1
!                    ELSE
                      EXIT TOPCELL
!                    END IF
                  END DO TOPCELL
                END IF
                IF ( ilay.LT.nlayp1 ) THEN
                  il = ilay
                  h = HNEW(ic, ir, il)
                ELSE
                  h = DBLE(BOTM(ic,ir,NLAY))
                END IF
                land = ABS(IUZFBND(ic, ir))
crgn changed HNEW(ic, ir, il) to h in next line.
                HLDUZF(ic, ir) = h
                IF ( IBOUND(ic, ir, il).LT.0 ) IUZFBND(ic, ir) = 0
                IF ( IUZFOPT.GT.0 ) THEN
C
C17-----SET CELL TOP, LENGTH, WIDTH AND WATER TABLE ELEVATION.
                  slen = DELC(ir) 
                  width = DELR(ic) 
! RGN changed BOTM(ic, ir, 0) to BOTM(ic, ir, land-1) 1/28/2010
                  IF ( land.GT.0 ) THEN
                    celtop = BOTM(ic, ir, land-1) - 0.5 * SURFDEP
                  ELSE
                    celtop = BOTM(ic, ir, 0) - 0.5 * SURFDEP
                  END IF
C
C18-----SKIP IF CELL IS OUTSIDE ACTIVE BOUNDARY OR IS NOT WATER TABLE.
! commented next line out to simulate unsat. flow over a portio of area.
!                  IF ( il.LT.1 ) IUZFBND(ic, ir) = 0
C
C19-----INITIALIZE UZTHST ARRAY TO RESIDUAL WATER CONTENT.

                  thtrcell = THTR(ic, ir)
                  DO jk = 1, NWAV
                    UZTHST(jk, l) = thtrcell
                  END DO
C
C20-----INITIALIZE UNSATURATED ZONE ARRAYS FOR FIRST STRESS PERIOD.
                  IF ( iflginit.EQ.1 ) THEN
                    IF ( celtop.GT.h ) THEN
                      UZDPST(1, l) = (celtop-h)
C
C21-----CALCULATE INITIAL WATER CONTENT AND FLUX IF STEADY STATE.
                      IF ( iss.NE.0 ) THEN
                        IF ( ITHTIFLG.GT.0 ) THEN
                          top = THTI(ic, ir) - thtrcell
                          IF ( top.LT.0.0 ) top = 0.0
                          UZFLST(1, l) = VKS(ic,ir)*
     +                       (top/(THTS(ic,ir)-thtrcell))**EPS(ic,ir)
                          UZTHST(1, l) = THTI(ic, ir)
                        ELSE
                          UZFLST(1, l) = FINF(ic, ir)
                          UZTHST(1, l) = (((UZFLST(1, l)/VKS(ic,ir))**
     +                    (1.0/EPS(ic,ir)))*(THTS(ic,ir)-thtrcell))
     +                                 + thtrcell
                          top = UZTHST(1, l) - thtrcell  
                        END IF
                        IF ( UZTHST(1, l)-thtrcell.LT.0.0D0 )
     +                       UZTHST(1, l) = thtrcell
C
C22-----SET INITIAL WATER CONTENT TO THTI AND CALCULATE FLUX IF 
C         TRANSIENT.
                      ELSE
                        UZTHST(1, l) = THTI(ic, ir)
                        top = UZTHST(1, l) - thtrcell
                        IF ( top.LT.0.0 ) top = 0.0
                        IF ( top.GT.0.0 ) THEN
                          bottom = THTS(ic, ir) - thtrcell
                          UZFLST(1, l) = VKS(ic, ir)*(top/bottom)
     +                                   **EPS(ic, ir)
                        END IF
                      END IF
                      IF ( UZTHST(1, l).LT.thtrcell ) UZTHST(1, l)
     +                     = thtrcell
C
C23-----CALCULATE VOLUME OF WATER STORED IN UNSATURATED ZONE.
                      IF ( top.GT.0.0 ) THEN
                        IF ( iss.EQ.0 ) UZSTOR(ic, ir) = UZDPST(1, l)
     +                       *top*width*slen
                        UZSPST(1, l) = 0.0D0
                        UZOLSFLX(ic, ir) = UZFLST(1, l)
                      ELSE
                        UZSTOR(ic, ir) = 0.0D0
                        UZFLST(1, l) = 0.0D0
                        UZSPST(1, l) = 0.0D0
                        UZOLSFLX(ic, ir) = 0.0D0
                      END IF
C
C24-----IF NO UNSATURATED ZONE, SET ARRAY VALUES TO ZERO EXEPT WHEN
C         STEADY STATE, THEN SET UZFLST ARRAY TO INFILRATION RATE.
                    ELSE
                      IF ( iss.NE.0 ) THEN
                        UZFLST(1, l) = FINF(ic, ir)
                      ELSE
                        UZFLST(1, l) = 0.0D0
                      END IF
                      UZDPST(1, l) = 0.0D0
                      UZSPST(1, l) = 0.0D0
                      UZTHST(1, l) = thtrcell
                      UZSTOR(ic, ir) = 0.0D0
cupdate        
                      UZOLSFLX(ic, ir) = FINF(ic, ir)
                    END IF
                    IF( RTSOLUTE.GT.0 ) THEN
                      DO uzlay = 1, NLAY
                        IF ( IBOUND(ic, ir, uzlay).LT.1 .OR.
     +                      HNEW(ic, ir, uzlay).LT.
     +                      BOTM(ic, ir, uzlay)) THEN
                           GRIDSTOR(ic, ir, uzlay) = 
     +                          (UZTHST(1, l)-thtrcell)*
     +                          (BOTM(ic,ir,uzlay-1)-BOTM(ic,ir,uzlay))
                        ELSEIF ( HNEW(ic, ir, uzlay).LT.
     +                           BOTM(ic, ir, uzlay-1)) THEN
                           GRIDSTOR(ic, ir, uzlay) = 
     +                        (UZTHST(1, l)-thtrcell)*
     +                        (BOTM(ic,ir,uzlay-1)-HNEW(ic, ir, uzlay))
                        ELSE
                          GRIDSTOR(ic, ir, uzlay) = 
     +                        (THTS(IC, IR)-thtrcell)*
     +                        (BOTM(ic,ir,uzlay-1)-BOTM(ic,ir,uzlay))
                        END IF
                      END DO
                    END IF
C
C25-----INITIALIZE ARRAYS FOR A TRANSIENT PERIOD THAT FOLLOWS A
C         STEADY STATE PERIOD IN STRESS PERIOD 1.
                  ELSE IF ( iflginit.EQ.2 ) THEN
                    IF ( celtop.GT.h ) THEN
                      UZDPST(1, l) = celtop - h
                      HLDUZF(ic,ir) = h
C
C26-----CALCULATE INITIAL WATER CONTENT AND FLUX FROM STEADY STATE
C         SIMULATION.
                      IF ( UZFLST(1, l).LT.0.0D0 ) UZFLST(1, l) = 0.0D0
                      IF ( ITHTIFLG.GT.0 ) THEN
                        top = THTI(ic, ir) - thtrcell
                        IF ( top.LT.0.0 ) top = 0.0
                        UZFLST(1, l) = VKS(ic,ir)*
     +                     (top/(THTS(ic,ir)-thtrcell))**EPS(ic,ir)
                        UZTHST(1, l) = THTI(ic, ir)
                      ELSE
                        UZTHST(1, l) = (((UZFLST(1, l)/VKS(ic,ir))**
     +                    (1.0/EPS(ic,ir)))*(THTS(ic,ir)-thtrcell))
     +                             + thtrcell
                      END IF
                      IF ( UZTHST(1, l).LT.thtrcell ) UZTHST(1, l)
     +                     = thtrcell
                      top = UZTHST(1, l) - thtrcell
                      IF ( top.LT.0.0 ) top = 0.0
                      IF ( top.LT.1.0E-5 ) UZFLST(1, l) = 0.0D0
                      IF ( top.GT.1.0E-5 ) THEN
                        UZSTOR(ic, ir) = UZDPST(1, l)*top*width*slen
                        UZSPST(1, l) = 0.0D0
                        UZOLSFLX(ic, ir) = UZFLST(1, l)
C
C27-----IF NO UNSATURATED ZONE, SET ARRAYS VALUES TO ZERO.
                      ELSE
                        UZSTOR(ic, ir) = 0.0D0
                        UZFLST(1, l) = 0.0D0
                        UZSPST(1, l) = 0.0D0
                        UZOLSFLX(ic, ir) = 0.0D0
                      END IF
                    ELSE
                      UZDPST(1, l) = 0.0D0
                      UZFLST(1, l) = 0.0D0
                      UZSPST(1, l) = 0.0D0
                      UZTHST(1, l) = thtrcell
                      UZSTOR(ic, ir) = 0.0D0
                      UZOLSFLX(ic, ir) = 0.0D0
                    END IF
                    IF( RTSOLUTE.GT.0 ) THEN
                      DO uzlay = 1, NLAY
                        IF ( IBOUND(ic, ir, uzlay).LT.1 .OR.
     +                      HNEW(ic, ir, uzlay).LT.
     +                      BOTM(ic, ir, uzlay)) THEN
                           GRIDSTOR(ic, ir, uzlay) = 
     +                          (UZTHST(1, l)-thtrcell)*
     +                          (BOTM(ic,ir,uzlay-1)-BOTM(ic,ir,uzlay))
                        ELSEIF ( HNEW(ic, ir, uzlay).LT.
     +                           BOTM(ic, ir, uzlay-1)) THEN
                           GRIDSTOR(ic, ir, uzlay) = 
     +                        (UZTHST(1, l)-thtrcell)*
     +                        (BOTM(ic,ir,uzlay-1)-HNEW(ic, ir, uzlay))
                        ELSE
                          GRIDSTOR(ic, ir, uzlay) = 
     +                        (THTS(IC, IR)-thtrcell)*
     +                        (BOTM(ic,ir,uzlay-1)-BOTM(ic,ir,uzlay))
                        END IF
                      END DO
                    END IF
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END IF
C
C28-----RETURN.
      RETURN
      END SUBROUTINE GWF2UZF1RP
!
!
      end module GwfUzfSubs

