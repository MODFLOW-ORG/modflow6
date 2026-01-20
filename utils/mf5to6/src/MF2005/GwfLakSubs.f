      module GwfLakSubs
        
        use GWFLAKMODULE, only: SGWF2LAK7PNT
        use SimPHMFModule, only: store_warning, store_error
        private
        public :: GWF2LAK7AR, GWF2LAK7RP

      contains

      SUBROUTINE GWF2LAK7AR(IN,IUNITSFR,IUNITGWT,IUNITUZF,NSOL,IGRID)
C
C------USGS VERSION 7.1; JUNE 2006 GWF2LAK7AR;
C------UPDATED FOR MF-2005, FEBRUARY 6, 2012
C     ******************************************************************
C     INITIALIZE POINTER VARIABLES USED BY SFR1 TO SUPPORT LAKE3 AND
C     GAGE PACKAGES AND THE GWT PROCESS
C     ******************************************************************
C
      USE GWFLAKMODULE
      USE GLOBAL,       ONLY: IOUT, NCOL, NROW, NLAY, IFREFM, ITRSS,
     +                        NODES
      USE GWFSFRMODULE, ONLY: NSS
      use utl7module, only: URDCOM, URWORD, U2DINT, U2DREL 
C
C      ******************************************************************
C      ALLOCATE ARRAY STORAGE FOR LAKES
C      ******************************************************************
C
C      ------------------------------------------------------------------
C      SPECIFICATIONS:
       CHARACTER (LEN=40):: CARD
       CHARACTER*200 line
       double precision :: r
C      ------------------------------------------------------------------
Crsr  Allocate lake variables used by SFR even if lakes not active so that
C       argument lists are defined
      ALLOCATE (NLAKES, NLAKESAR,THETA,LAKUNIT)
      allocate (NeedLakWaterMover)
      NeedLakWaterMover = .false.
      NLAKES = 0
      LAKUNIT = IN
      NLAKESAR = 1
      THETA = 0.0
C0--If LAK package is active
      IF (IN.GT.0) THEN
Cdep added SURFDEPTH 3/3/2009
        ALLOCATE (ILKCB, NSSITR, SSCNCR, SURFDEPTH)
        ALLOCATE (MXLKND, LKNODE, ICMX, NCLS, LWRT, NDV, NTRB)
        ALLOCATE (IRDTAB)
C
C1------IDENTIFY PACKAGE AND INITIALIZE LKNODE.
      WRITE(IOUT,1) IN
      LKNODE=0
Cdep  initialize number of iterations and closure criteria to zero.
      DUM = 0.0
      NSSITR = 0
      SSCNCR = 0.0
      SURFDEPTH = 0.0
!
      lloc = 1
      IRDTAB = 0
      NPP = 0
      MXVL = 0
! Read item 1a
      CALL URDCOM(In, IOUT, line)
! Check for alternate option to specifiy stage/vol/area tables.
      CALL UPARLSTAL(IN,IOUT,LINE,NPP,MXVL)  ! ERB - Pointless, since LAK does not support parameters
      lloc = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'TABLEINPUT') THEN
         IRDTAB = 1
         WRITE(IOUT,32)
   32  FORMAT(1X,I10,' Stage, volume and area relationship specified ',
     +                'based on an external tabular input file')
      ELSE
        BACKSPACE IN
        WRITE(IOUT,'(A)') ' Model grid will be used to develop ',
     +                     ' volume and area relationship. '
      ENDIF
C
C2------READ NLAKES, ILKCB.
C
Cdep  Revised input statement to read THETA,NSSITR,SSCNCR for
Cdep  transient simulations when THETA is negative.
        IF(IFREFM.EQ.0) THEN
! Read item 1b
           READ(IN,'(2I10)')NLAKES,ILKCB
! Read item 2 and backspace
           IF (ITRSS.LE.0) THEN
              READ(IN,'(F10.2,I10,F10.2)') THETA,NSSITR,SSCNCR
              IF (THETA.LT.0.0) BACKSPACE IN
           ELSE
              READ(IN,'(F10.2)') THETA
              IF (THETA.LT.0.0) BACKSPACE IN
           ENDIF
        ELSE
! Read item 1b
           READ(IN,*) NLAKES,ILKCB
! Read item 2 and backspace
           IF (ITRSS.LE.0) THEN
              READ(IN,*) THETA,NSSITR,SSCNCR
              IF(THETA.LT.0.0) BACKSPACE IN
           ELSE
              READ(IN,*) THETA
              IF(THETA.LT.0.0) BACKSPACE IN
           ENDIF
        ENDIF

Cdep    Set default values for number of iterations and closure criteria
Cdep     for transient simulations when using original version of
Cdep     LAKE Package.
        IF(THETA.GE.0.0.AND.NSSITR.EQ.0) THEN
          NSSITR=100
          SSCNCR=1.0E-05
        ELSEIF(THETA.LT.0.0)THEN
          THETA=ABS(THETA)
! Read item 2
          IF(IFREFM.EQ.0) THEN
Cdep fixed format can't read in exponent notation
!rsr, old data sets may not have SURFDEPTH, may need to trap this for some compilers
            READ (IN, '(A)') CARD
            NUMCHAR = LEN(TRIM(CARD))
            IF ( NUMCHAR>30 ) THEN
              READ(CARD,'(F10.2,I10,2F10.5)') DUM,NSSITR,SSCNCR,
     +                                         SURFDEPTH
            ELSE
              READ(CARD,'(F10.2,I10,F10.5)') DUM,NSSITR,SSCNCR
            ENDIF
          ELSE
            READ(IN,*,IOSTAT=IOS) DUM,NSSITR,SSCNCR,SURFDEPTH
            IF ( IOS.NE.0 ) SURFDEPTH = 0.0
          ENDIF
        ENDIF
Cdep   Add check to reset THETA when > 1 or < 0.5.
        IF(THETA.GT.1.0) THEN
          THETA = 1.0
        ELSEIF(THETA.LT.0.5)THEN
          THETA = 0.0
        ENDIF
      ENDIF   ! goes with IF at comment C0
C
C
C  SET NLAKES ARRAY VARIABLE TO NLAKES IF NLAKES GREATER THAN 0.
      IF (NLAKES.GT.0) NLAKESAR = NLAKES
      ALLOCATE (VOL(NLAKESAR), STGOLD(NLAKESAR), STGNEW(NLAKESAR))
      ALLOCATE(STGOLD2(NLAKESAR))
      ALLOCATE (VOLOLDD(NLAKESAR))
!     ALLOCATE (VOLOLDD(NLAKESAR), VOLOLD(NLAKES), VOLINIT(NLAKES))
      ALLOCATE (STGITER(NLAKESAR))
      ALLOCATE (LAKSEEP(NCOL,NROW))
      STGNEW = 0.0D0
      STGOLD = 0.0D0
      STGOLD2 = 0.0D0
      STGITER = 0.0D0
      VOLOLDD = 0.0D0
      LAKSEEP = 0.0
Cdep initialized VOLOLD and VOLINIT  6/4/2009 (VOLOLD is single precision)
!     VOLOLD = 0.0
!     VOLINIT = 0.0
      VOL = 0.0
      CALL SGWF2LAK7PSV1(IGRID)
      IF (IN.LT.1) RETURN
C
C Lakes are active
      ALLOCATE (STAGES(NLAKESAR), CLAKE(NLAKESAR,NSOL))
      STAGES = 0.0
      CLAKE = 0.0
C Budget variables for GSFLOW
      ALLOCATE (TOTGWIN_LAK,TOTGWOT_LAK,TOTDELSTOR_LAK,TOTSTOR_LAK)
      ALLOCATE (TOTEVAP_LAK,TOTPPT_LAK,TOTRUNF_LAK,TOTWTHDRW_LAK)
      ALLOCATE (TOTSURFIN_LAK,TOTSURFOT_LAK)
      TOTGWIN_LAK = 0.0
      TOTGWOT_LAK = 0.0
      TOTDELSTOR_LAK = 0.0
      TOTSTOR_LAK = 0.0
      TOTEVAP_LAK = 0.0
      TOTPPT_LAK = 0.0
      TOTRUNF_LAK = 0.0
      TOTWTHDRW_LAK = 0.0
      TOTSURFIN_LAK = 0.0
      TOTSURFOT_LAK = 0.0
C
C  VALUE OF MXLKND (NUMBER OF LAKE-AQUIFER INTERFACES) IS AN ESTIMATE.
C    TO SAVE MEMORY, REDUCE ITS SIZE IF APPROPRIATE.
C    IF MXLKND TOO SMALL, ERROR MESSAGE WILL BE PRINTED.
      MXLKND=NCOL*NROW*NLAY/2
      IF (NLAKES.LT.1) THEN
        WRITE(IOUT,2)
        IN=0
        NLAKES = 0
      ELSE
      WRITE(IOUT,5) MXLKND,NLAKES
      IF (ILKCB.GT.0) WRITE(IOUT,7) ILKCB
      IF (ILKCB.LE.0) WRITE(IOUT,9)
Cdep   Write THETA, NSSITR, SSCNCR
      IF (ITRSS.GT.0) THEN
        WRITE(IOUT,22) THETA
        WRITE(IOUT,10) NSSITR, SSCNCR
      ELSE
        WRITE(IOUT,11) THETA, NSSITR, SSCNCR
      ENDIF
Cdep   Changed default values for NSSITR and SSCNCR and revised
Cdep     print statements using format statement 10.
Cdep      IF(ITRSS.LE.0.AND.NSSITR.EQ.0) NSSITR = 50
Cdep      IF(ITRSS.LE.0.AND.SSCNCR.EQ.0.0) SSCNCR = 0.01
Cdep      IF(ITRSS.EQ.0) WRITE(IOUT,23) NSSITR, SSCNCR
Cdep      IF(ITRSS.LT.0) WRITE(IOUT,24) NSSITR, SSCNCR
C-lfk  1     FORMAT(/1X,'LAK7 -- LAKE PACKAGE, VERSION 7, 2/06/2012',
1     FORMAT(/1X,'LAK7 -- LAKE PACKAGE, VERSION 7, 1/07/2013',
     1' INPUT READ FROM UNIT',I3)
2       FORMAT(1X,' NUMBER OF LAKES=0, ',
     1              ' SO LAKE PACKAGE IS BEING TURNED OFF')
5     FORMAT(1X,'SPACE ALLOCATION FOR',I7,' GRID CELL FACES ADJACENT TO
     1LAKES'/1X,'MAXIMUM NUMBER OF LAKES IS',I3, ' FOR THIS SIMULATION')
7     FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE RECORDED ON UNIT',I5)
9     FORMAT(1X,'CELL-BY-CELL SEEPAGES WILL NOT BE PRINTED OR SAVED')
Cdep added format statement when starting with transient simulation
  10  FORMAT(//1X,'LAKE PACKAGE HAS BEEN MODIFIED TO ITERATIVELY ',
     1 'SOLVE FOR LAKE STAGE DURING TRANSIENT STRESS PERIODS:',/1X,
     2 'MAXIMUM NUMBER OF ITERATIONS (NSSITR) = ',I5,/1X,
     3 'CLOSURE CRITERIA FOR LAKE STAGE (SSCNCR) = ',1PE12.6,/1X,
     4 'DEFAULT VALUES FOR TRANSIENT ONLY SIMULATIONS ARE: ',
     5 'NSSITR = 100 AND SSCNCR = 0.0001',/1X,'VALUES OTHER THAN ',
     6 'DEFAULT CAN BE READ BY SPECIFYING A THETA LESS THAN ZERO ',
     7 'THEN ADDING NSSITR AND SSCNCR PER ORIGINAL INSTRUCTIONS.',/1X,
     8 'NEGATIVE THETA MUST BE LESS THAN ZERO BUT NOT MORE THAN ',
     9 'ONE. THETA IS CONVERTED TO A POSITIVE VALUE.',/1X,
     * 'MINIMUM AND MAXIMUM LAKE STAGES FOR TRANSIENT ',
     * 'SIMULATIONS ARE SET TO BOTTOM AND TOP ELEVATIONS USED TO ',
     * 'COMPUTE LAKE VOLUME, RESPECTIVELY.',//)
Cdep added format statement for steady state only simulations.
  11  FORMAT(//1X,'NEWTON ITERATION METHOD FOR COMPUTING LAKE STAGE ',
     1 'DURING STEADY-STATE STRESS PERIODS HAS BEEN MODIFIED:',/1X,
     2 'SPECIFIED THETA OF ',F6.3,' WILL BE AUTOMATICALLY CHANGED TO ',
     3 '1.0 FOR ALL STEADY STATE STRESS PERIODS.',/1X,
     4 'MAXIMUM NUMBER OF STEADY-STATE ITERATIONS (NSSITR) = ',I5,/1X,
     5 'CLOSURE CRITERIA FOR STEADY-STATE LAKE STAGE (SSCNCR) = ',
     6  1PE12.6,//)
Cdep revised print statement to note that time weighting of theta can
Cdep  vary only between 0.5 and 1 for transient simulations
Cdep   22 FORMAT(/1X,'THETA = ',F10.2,'  METHOD FOR UPDATING LAKE STAGES IN
Cdep     1ITERATIONS OF THE SOLUTION FOR AQUIFER HEADS.'/20X,'0.0 IS EXPLICI
Cdep     2T, 0.5 IS CENTERED, AND 1.0 IS FULLY IMPLICIT.')
   22 FORMAT(/1X,'THETA = ',F6.3,/1X,'THETA IS THE TIME WEIGHTING ',
     *'FACTOR FOR COMPUTING LAKE STAGE DURING TRANSIENT MODFLOW ',
     *'TIME STEPS AND ITS DEFINITION HAS BEEN MODIFIED.',/1X,'A THETA ',
     *'OF LESS THEN 0.5 IS AUTOMATICALLY SET TO 0 AND LAKE STAGE IS ',
     *'EQUAL TO THE STAGE AT THE END OF THE PREVIOUS TIME STEP. ',/1X,
     *'TRANSIENT SIMULATIONS OF LAKE STAGE WITH THE CURRENT TIME STEP ',
     *'REQUIRES A THETA BETWEEN 0.5 AND 1.0. ',/1X,'VALUES GREATER ',
     *'THAN 1.0 ARE AUTOMATICALLY RESET TO  1.0 AND VALUES LESS ',
     *'THAN 0.5 ARE RESET TO 0.0.',/1X,'A THETA OF 0.5 REPRESENTS THE ',
     *'AVERAGE LAKE STAGE DURING A TIME STEP.',/1X,'A THETA OF 1.0 ',
     *'REPRESENTS THE LAKE STAGE AT THE END OF THE TIME STEP.',//)
Cdep   23 FORMAT(/1X,'STEADY-STATE SOLUTION FOR LAKES.'
Cdep     2/1X,'MAXIMUM NUMBER OF ITERATIONS = ',I4,3X,
Cdep     1'CONVERGENCE CRITERION = ',1PE9.2)
Cdep   24 FORMAT(/1X,'COMBINED STEADY-STATE/TRANSIENT SOLUTION FOR LAKES.'
Cdep     2/1X,'MAXIMUM NUMBER OF ITERATIONS = ',I4,3X,
Cdep     1'CONVERGENCE CRITERION = ',1PE9.2)

        ALLOCATE (ILAKE(5,MXLKND), BEDLAK(MXLKND), CNDFCT(MXLKND))
        ALLOCATE (PRCPLK(NLAKES), EVAPLK(NLAKES), WTHDRW(NLAKES))
        ALLOCATE (RNF(NLAKES), CRNF(NLAKES,NSOL), CUMRNF(NLAKES))
        ALLOCATE (CUMUZF(NLAKES))
        ALLOCATE (ISUB(NLAKES,NLAKES), SILLVT(NLAKES,NLAKES))
        ALLOCATE (IRK(2,NLAKES))
        ALLOCATE (CUMPPT(NLAKES), CUMEVP(NLAKES), CUMGWI(NLAKES))
        ALLOCATE (CUMGWO(NLAKES), CUMSWI(NLAKES), CUMSWO(NLAKES))
        ALLOCATE (CUMWDR(NLAKES), CUMFLX(NLAKES))
        ALLOCATE (CAUG(NLAKES,NSOL), CPPT(NLAKES,NSOL))
        ALLOCATE (CLAKINIT(NLAKESAR,NSOL))
        ALLOCATE (ICS(NLAKES),BOTTMS(NLAKES), BGAREA(NLAKES))
        ALLOCATE (SSMN(NLAKES), SSMX(NLAKES))
        ALLOCATE (LKARR1(NCOL,NROW,NLAY), BDLKN1(NCOL,NROW,NLAY))
        ALLOCATE (EVAP(NLAKES), PRECIP(NLAKES), SEEP(NLAKES),
     +            SEEP3(NLAKES),EVAP3(NLAKES), PRECIP3(NLAKES))
        ALLOCATE (SEEPUZ(NLAKES))
        ALLOCATE (FLWITER(NLAKES),FLWITER3(NLAKES))
        ALLOCATE (SURFA(NLAKES), SURFIN(NLAKES), SURFOT(NLAKES))
        ALLOCATE (SUMCNN(NLAKES), SUMCHN(NLAKES))
        ALLOCATE (NCNCVR(NLAKES), LIMERR(NLAKES), DSRFOT(NLAKES))
Cdep  Allocate arrays that track lake budgets for dry lakes
        ALLOCATE (EVAPO(NLAKES),WITHDRW(NLAKES),FLWIN(NLAKES))
        ALLOCATE (GWRATELIM(NLAKES))
        EVAPO = 0.0
        WITHDRW = 0.0D0
        FLWIN = 0.0
        FLWITER = 0.0D0
        FLWITER3 = 0.0D0
        EVAP = 0.0D0
        PRECIP = 0.0D0
        EVAP3 = 0.0D0
        PRECIP3 = 0.0D0
        IF ( IRDTAB.GT.0 ) THEN
          ALLOCATE(LAKTAB(NLAKES))
        ELSE
          ALLOCATE(LAKTAB(1))
        ENDIF
        LAKTAB = 0
!rsr    GWRATLIM= 0.0
Cdep  Allocate space for three arrays used in GAGE Package
C       when Solute Transport is active
        ALLOCATE (XLAKES(NLAKES,1), XLAKINIT(NLAKES,1))
        ALLOCATE (XLKOLD(NLAKES,1))
crsr  Allocate arrays for BD subroutine
        ALLOCATE (LDRY(NODES), FLXINL(NLAKES))
        ALLOCATE (NCNT(NLAKES), NCNST(NLAKES))
        ALLOCATE (SVT(NLAKES), KSUB(NLAKES), STGADJ(NLAKES))
        ALLOCATE (MSUB(NLAKES,NLAKES), MSUB1(NLAKES))
        ALLOCATE (GWIN(NLAKES), GWOUT(NLAKES))
        ALLOCATE (DELH(NLAKES), TDELH(NLAKES))
Cdep   Allocate lake budget error arrays for BD subroutine 6/9/2009
        ALLOCATE (CUMVOL(NLAKES), CMLAKERR(NLAKES))
        ALLOCATE (CUMLKIN(NLAKES), CUMLKOUT(NLAKES))
        ALLOCATE (DELVOL(NLAKES), TSLAKERR(NLAKES))
Cdep initialized VOLOLD and VOLINIT  6/4/2009 (VOLOLD is single precision)
        ALLOCATE (VOLOLD(NLAKES), VOLINIT(NLAKES))
        VOLOLD = 0.0
        VOLINIT = 0.0
      ENDIF
Cdep   ALLOCATE SPACE FOR CONNECTION WITH STREAMS
      IF (IUNITSFR.LE.0) THEN
        NSSAR = 1
      ELSE
        NSSAR = NSS
      ENDIF
Cdep   ALLOCATE SPACE FOR FLOB ARRAY WHEN TRANSPORT ACTIVE.
      IF (IUNITGWT.LE.0) THEN
        MXLKAR = 1
      ELSE
        MXLKAR = MXLKND
      ENDIF
Cdep    ALLOCATE SPACE FOR OVERLAND FLOW WHEN UNSATURATED FLOW ACTIVE.
! RGN Allocate NUZFAR to nlakes for all cases because of the GAG package 5/28/09
!      IF (IUNITUZF.LE.0) THEN
!       NUZFAR = 1
!      ELSE
        NUZFAR = NLAKESAR
!      ENDIF

      !rsr, what if NLAKES < 1, sanity check
      IF (NLAKES<1 ) THEN
        print *, 'nlakes dimension problem in lak7', nlakes
        stop
      ENDIF

      ALLOCATE (ITRB(NLAKES,NSSAR), IDIV(NLAKES,NSSAR))
      ALLOCATE (FLOB(MXLKAR))
      ALLOCATE (OVRLNDRNF(NUZFAR), CUMLNDRNF(NUZFAR))
Cdep    ALLOCATE SPACE FOR DEPTHTABLE, AREATABLE, AND VOLUMETABLE
      ALLOCATE (DEPTHTABLE(151,NLAKES), AREATABLE(151,NLAKES))
      ALLOCATE (VOLUMETABLE(151,NLAKES))
      ITRB = 0
      IDIV = 0
      FLOB = 0.0
      OVRLNDRNF = 0.0
      CUMLNDRNF = 0.0
      CUMUZF = 0.0
      DEPTHTABLE = 0.0D0
      AREATABLE = 0.0D0
      VOLUMETABLE = 0.0D0
Cdep initialized lake budget error arrays  6/9/2009
      CUMVOL = 0.0
      DELVOL = 0.0
      CMLAKERR = 0.0
      TSLAKERR = 0.0
      CUMLKOUT = 0.0
      CUMLKIN = 0.0
C-----SAVE POINTERS FOR GRID AND RETURN
      CALL SGWF2LAK7PSV(IGRID)
C
C11-----RETURN.
      RETURN
      END SUBROUTINE GWF2LAK7AR
C
      SUBROUTINE GWF2LAK7RP(IN,IUNITBCF,IUNITGWT,IUNITLPF,IUNITHUF,
     +                      IUNITSFR,IUNITUZF,IUNITUPW,KKPER,NSOL,
     +                      IOUTS,IGRID)
C
C------USGS VERSION 7.1;  JUNE 2006 GWF2LAK7RP
C        REVISED FEBRUARY 6, 2012
C     ******************************************************************
C       READ INPUT DATA FOR THE LAKE PACKAGE.
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLAKMODULE
      USE GLOBAL,       ONLY: IOUT, NCOL, NROW, NLAY, IFREFM, IBOUND,
     +                        LBOTM, BOTM, DELR, DELC, ISSFLG
      use SimPHMFModule, only: ustop
      use utl7module, only: U2DINT, U2DREL
C     USE GWFSFRMODULE, ONLY: NSS
C     ------------------------------------------------------------------
C     FUNCTIONS
C     ------------------------------------------------------------------
!*      DOUBLE PRECISION VOLTERP
!*      EXTERNAL VOLTERP
C     ------------------------------------------------------------------
      CHARACTER*24 ANAME(2)
!     CHARACTER*30 LFRMAT
!dep  added STGINIT as double precision
      DOUBLE PRECISION STGINIT
      logical, save :: warned = .false.
      character(len=200) :: warning
      DATA ANAME(1)/'           LAKE ID ARRAY'/
      DATA ANAME(2)/'  LAKEBED LEAKANCE ARRAY'/
C
C     ------------------------------------------------------------------
C------SET POINTERS FOR THE CURRENT GRID.
      CALL SGWF2LAK7PNT(IGRID)
C
C1A-----IF MXLKND IS LESS THAN 1, THEN LAKE IS INACTIVE. RETURN.
      IF(MXLKND.LT.1) RETURN
C
C1A1----READ INITIAL CONDITIONS FOR ALL LAKES (ONLY READ ONCE)
      ISS = ISSFLG(KKPER)
      IF (KKPER.EQ.1) THEN
         WRITE (IOUT,19)
         IF(ISS.NE.0) WRITE (IOUT,20)
         IF(ISS.EQ.0) WRITE (IOUT,820)
! Read Item 3         
         IF (IUNITGWT.EQ.0) THEN
            DO 30 LM=1,NLAKES
               IF (IFREFM.EQ.0) THEN
                 IF ( IRDTAB.GT.0 ) THEN
                   IF(ISS.NE.0) READ (IN,'(3F10.4,I5)') STAGES(LM),
     1                                SSMN(LM),SSMX(LM),LAKTAB(LM)
                   IF(ISS.EQ.0) READ (IN,'(F10.4,I5)') STAGES(LM),
     2                                               LAKTAB(LM)
                 ELSE
                   IF(ISS.NE.0) READ (IN,'(3F10.4)') STAGES(LM),
     1              SSMN(LM),SSMX(LM)
                   IF(ISS.EQ.0) READ (IN,'(F10.4)') STAGES(LM)
                 ENDIF
               ELSE
                 IF ( IRDTAB.GT.0 ) THEN
                   IF(ISS.NE.0) READ (IN,*)STAGES(LM),SSMN(LM),SSMX(LM),
     1                                     LAKTAB(LM)
                   IF(ISS.EQ.0) READ (IN,*) STAGES(LM),LAKTAB(LM)
                 ELSE
                   IF(ISS.NE.0) READ (IN,*) STAGES(LM),SSMN(LM),SSMX(LM)
                   IF(ISS.EQ.0) READ (IN,*) STAGES(LM)
                 ENDIF
               ENDIF
            IF(ISS.NE.0) WRITE (IOUT,22) LM,STAGES(LM),SSMN(LM),SSMX(LM)
            IF(ISS.EQ.0) WRITE (IOUT,22) LM,STAGES(LM)
 30         CONTINUE
         ELSE
            WRITE (IOUTS,21) NSOL
!            WRITE (LFRMAT,23) NSOL  !LFRMAT is not set
            DO 35 LM=1,NLAKES
               IF (IFREFM.EQ.0) THEN
                 IF ( IRDTAB.GT.0 ) THEN
                   IF(ISS.NE.0) READ(IN,'(100F10.4)') STAGES(LM),
     1                SSMN(LM),SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL),
     2                LAKTAB(LM)
                   IF(ISS.EQ.0) READ (IN,'(100F10.4)') STAGES(LM),
     1                          (CLAKE(LM,ISOL),ISOL=1,NSOL),LAKTAB(LM)
                 ELSE
                   IF(ISS.NE.0) READ(IN,'(100F10.4)') STAGES(LM),
     1                SSMN(LM),SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL)
                   IF(ISS.EQ.0) READ (IN,'(100F10.4)') STAGES(LM),
     1                          (CLAKE(LM,ISOL),ISOL=1,NSOL)
                 ENDIF
               ELSE
                 IF ( IRDTAB.GT.0 ) THEN
                   IF(ISS.NE.0) READ (IN,*) STAGES(LM),SSMN(LM),
     1                          SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL),
     2                          LAKTAB(LM)
                   IF(ISS.EQ.0) READ (IN,*) STAGES(LM),
     1                          (CLAKE(LM,ISOL),ISOL=1,NSOL),LAKTAB(LM)
                 ELSE
                   IF(ISS.NE.0) READ (IN,*) STAGES(LM),SSMN(LM),
     1                          SSMX(LM),(CLAKE(LM,ISOL),ISOL=1,NSOL)
                   IF(ISS.EQ.0) READ (IN,*) STAGES(LM),
     1                          (CLAKE(LM,ISOL),ISOL=1,NSOL)
                 ENDIF
               ENDIF
            IF(ISS.NE.0) WRITE (IOUT,22) LM,STAGES(LM),SSMN(LM),SSMX(LM)
            IF(ISS.EQ.0) WRITE (IOUT,22) LM,STAGES(LM)
 35           WRITE (IOUTS,*) LM,(CLAKE(LM,ISOL),ISOL=1,NSOL)
cgage
C            CLAKINIT=CLAKE
         ENDIF
      ENDIF
C
      WRITE (IOUT,'(/)')
      WRITE(IOUT,822)
 19   FORMAT(//1X,'LAKE PACKAGE ACTIVE:  CALCULATED LAKE STAGE FOR EACH
     1TIME STEP WILL BE STORED IN HNEW ARRAY.')
 20   FORMAT(///1X,'INITIAL LAKE STAGE:  LAKE    STAGE    SS MIN    SS M
     1AX'/)
 21   FORMAT (//1X,'INITIAL LAKE CONCENTRATIONS:  LAKE   CONCENTRATION (
     1NSOL =',I3,')'/)
 22   FORMAT (22X,I3,3F10.3)
 23   FORMAT ('(31X,I3,3X,1P',I3,'(E12.3))')
 820  FORMAT (/1X,'INITIAL LAKE STAGE:  LAKE    STAGE'/)
 822  FORMAT(//1X,'If any subsequent steady-state stress periods, min. a
     1nd max. stages for each lake will be read in Record 9a.'//)
C
! RGN 9/25/12 moved this to read lake bathymetry before stress period information.
      IF ( KKPER==1 .AND. IRDTAB.GT.0 ) THEN
        DO L1=1,NLAKES
          WRITE(IOUT,1399) L1
          iunit = LAKTAB(L1)
 1399 FORMAT(//1X,'STAGE/VOLUME RELATION FOR LAKE',I3//6X,'STAGE',
     1        8X,'VOLUME',8X,'AREA'/)
          DO  INC=1,151
          READ(iunit,*) DEPTHTABLE(INC,L1), VOLUMETABLE(INC,L1),
     +                    AREATABLE(INC,L1)
          WRITE(IOUT,1315) DEPTHTABLE(INC,L1), VOLUMETABLE(INC,L1),
     +                    AREATABLE(INC,L1)
          ENDDO
        ENDDO
      ENDIF
C1B-----READ ITMP (FLAG TO REUSE LAKE-GEOMETRY DATA).
! Read Item 4
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(3I10)') ITMP, ITMP1, LWRT
      ELSE
         READ(IN,*) ITMP, ITMP1, LWRT
      ENDIF
C
C2A-----IF ITMP < 0 THEN REUSE LAKE CONFIGURATION DATA FROM LAST STRESS
C       PERIOD.
      IF(ITMP.GE.0) GO TO 50
      WRITE (IOUT,'(/)')
      WRITE(IOUT,2)
    2 FORMAT(1H ,'REUSING LAKE CONFIGURATION DATA FROM LAST STRESS PERIO
     1D'/)
      GO TO 800
C
C4------IF THERE ARE NO LAKE NODES THEN RETURN.
   50 LKNODE = 0
      IF(ITMP.EQ.0) GOTO 900
      if (KKPER > 1 .and. .not. warned) then
        ! write warning about changing lakes
   55 format('In LAK input, ITMP > 0 for stress period ',i0,
     & '.  WARNING: LAK8 does not support changing lake configuration',
     & ' or leakance during simulation.')
        write(warning,55)kkper
        call store_warning(warning)
        warned = .true.
      endif
C
C   INITIALIZE BGAREA
      DO 60 LK=1,NLAKES
      BGAREA(LK)=0.0
   60 CONTINUE
C
C5------READ INTEGER ARRAYS THAT DEFINE THE POSITIONS OF ALL LAKES IN
C5A     EACH MODEL GRID LAYER.  THEN READ ARRAYS OF LAKEBED CONDUCTANCES
C5B     IN EACH LAYER.
C
C   READ ARRAY OF LAKE ID'S, LAYER BY LAYER
C   REVISED 11/30/2005 DEP
! Read item 5
      DO 125 K=1,NLAY
      KK = K
      CALL U2DINT(LKARR1(:,:,KK),ANAME(1),NROW,NCOL,KK,IN,IOUT)
  125 CONTINUE
C
C   CHECK THAT ALL ENTRIES ARE VALID LAKE ID NUMBERS OR ZERO
C
      DO 130 K=1,NLAY
      DO 130 I=1,NCOL
      DO 130 J=1,NROW
      IF(LKARR1(I,J,K).GT.0.AND.LKARR1(I,J,K).LE.NLAKES) GO TO 130
      LKARR1(I,J,K)=0
  130 CONTINUE
C
C   CHECK IF LAKE CELLS HAVE VALUES OF IBOUND=0; WARN IF INCONSISTENT
C
      WRITE (IOUT,'(/)')
      DO 132 K=1,NLAY
      DO 132 I=1,NCOL
      DO 132 J=1,NROW
      IF(LKARR1(I,J,K).GT.0.AND.IBOUND(I,J,K).NE.0) THEN
         WRITE (IOUT,232) IBOUND(I,J,K),LKARR1(I,J,K),I,J,K
  232    FORMAT (7X,'*** WARNING: IBOUND = ',I2,
     1  ' & LKARR = ',I2,' at CELL I=',I3,
     2  ', J=',I3,', K=',I3,' ***')
      ENDIF
  132 CONTINUE
C
C   READ ARRAY OF BED LEAKANCES, LAYER BY LAYER
Cdep    REVISED 11/30/2005
      WRITE (IOUT,'(/)')
! Read item 6
      DO 135 K=1,NLAY
      KK = K
      CALL U2DREL(BDLKN1(:,:,KK),ANAME(2),NROW,NCOL,KK,IN,IOUT)
  135 CONTINUE
C
        WRITE(IOUT,36)
        WRITE(IOUT,4)
36    FORMAT(/7X,'LOCATIONS, LAKE #, INTERFACE TYPE FOR GRID CELLS',
     1 ' ADJACENT TO LAKES:',5X,/
     3 5X,71('-'))
4     FORMAT(5X,'LAYER #',4X,'ROW #',4X,'COLUMN #',3X,'LAKE #',
     1       2X,'INTERFACE TYPE',2X,'LAKEBED LEAKANCE')
C
C   IDENTIFY LAKE BORDER CELLS, ASSIGN CELL TYPE ID'S, COMPUTE AND
C     ASSIGN LAKE-AQUIFER INTERFACE CONDUCTANCES.
C
      M = 0
      DO 180 I=1,NCOL
      DO 180 J=1,NROW
      K = 1
      IF(LKARR1(I,J,K).EQ.0) GO TO 150
      IF(NLAY.EQ.1) GO TO 145
C   Keep searching in vertical direction until non-lake cell is found,
C     and define interface there ("K" for interface is layer below
C     bottom of lake)
      DO 140 K=2,NLAY
      IF(LKARR1(I,J,K).EQ.0) GO TO 145
  140 CONTINUE
C   Make sure that K=NLAY if lake extends to bottom cell of grid:
      K=NLAY
C      GO TO 145
C
C   VERTICAL LAKEBED INTERFACE (TYPE 0) DETECTED
C
  145 M = M + 1
      IF(M.LE.MXLKND) GO TO 147
      WRITE(IOUT,149) I,J,K
  149 FORMAT(/1X,'MAXIMUM NUMBER OF GRID CELLS ADJACENT TO LAKES HAS BEE
     1N EXCEEDED WITH CELL ',3I5,'  REDEFINE VARIABLE MXLKND TO A LARGER
     2 VALUE IN MODULE GWF2LAK7AR')
      CALL USTOP(' ')
  147 ILAKE(1,M) = K
      ILAKE(2,M) = J
      ILAKE(3,M) = I
Cdep  changed if statement August 24, 2009
Cdep      IF(K.GT.1.AND.LKARR1(I,J,K).EQ.0) LID = LKARR1(I,J,K-1)
Cdep      IF(LKARR1(I,J,K).NE.0) LID = LKARR1(I,J,K)
      IF(K.GT.1) THEN
        IF(LKARR1(I,J,K).EQ.0) THEN
          LID = LKARR1(I,J,K-1)
        ELSE
          LID = LKARR1(I,J,K)
        ENDIF
      ELSEIF (K.EQ.1) THEN
        IF(LKARR1(I,J,K).EQ.0) THEN
          LID = 0
        ELSE
          LID = LKARR1(I,J,K)
        ENDIF
      ENDIF
      ILAKE(4,M) = LID
      ILAKE(5,M) = 6
      IF ( K.GT.1 ) THEN             !RGN 5/21/12 added IF test
        BEDLAK(M) = BDLKN1(I,J,K-1)
      ELSE                           !RGN
        BEDLAK(M) = BDLKN1(I,J,K)    !RGN
      ENDIF                         !RGN
      IF(K.EQ.NLAY.AND.LKARR1(I,J,K).NE.0) BEDLAK(M) = 0.0
      BGAREA(LID) = BGAREA(LID) + DELC(J)*DELR(I)
C-LFK-JAN. 2013
c-lfk        WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
      WRITE(IOUT,6) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
5     FORMAT(5I10,10X,F10.5)
6     FORMAT(5I10,12X,1PE10.3)
C-LFK
      IF(LKARR1(I,J,K).NE.0) GO TO 180
C
C   SEARCH FOR CELL(S) ADJACENT TO LAKE
C
  150 K2 = K
      DO 175 K1=K2,NLAY
cgzh fix for 2D-problems
      IF(NCOL.EQ.1) GO TO 165
      IF(I.NE.1) GO TO 1151
      IF(LKARR1(I+1,J,K1).EQ.0) GO TO 165
      GO TO 1153
 1151 IF(I.NE.NCOL) GO TO 1152
      IF(LKARR1(I-1,J,K1).EQ.0) GO TO 165
      GO TO 1153
 1152 IF(LKARR1(I+1,J,K1).EQ.0.AND.LKARR1(I-1,J,K1).EQ.0) GO TO 165
C
C   CELL(S) LATERALLY ADJACENT TO LAKE IN X-DIRECTION (TYPE 1) DETECTED
C
 1153 DO 160 N=1,2
      IF(N.EQ.2) GO TO 155
      IF(I.EQ.1) GO TO 160
      IF(LKARR1(I-1,J,K1).EQ.0) GO TO 160
      I2 = I-1
      IFACE=1
      GO TO 157
  155 IF(I.EQ.NCOL) GO TO 160
      IF(LKARR1(I+1,J,K1).EQ.0) GO TO 160
      I2 = I + 1
      IFACE=2
  157 M = M + 1
      IF(M.LE.MXLKND) GO TO 158
      WRITE(IOUT,149) I,J,K1
      CALL USTOP(' ')
  158 ILAKE(1,M) = K1
      ILAKE(2,M) = J
      ILAKE(3,M) = I
      ILAKE(4,M) = LKARR1(I2,J,K1)
      ILAKE(5,M) = IFACE
      BEDLAK(M) = BDLKN1(I,J,K1)
      K4 = K1 - 1
      DO 3158 K3=1,K4
      IF(LKARR1(I,J,K3).EQ.0) GO TO 3158
      GO TO 3162
 3158 CONTINUE
      BEDLAK(M) = BDLKN1(I,J,1)
 3162 CONTINUE
c-lfk      WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
C-LFK-JAN. 2013
c-lfk        WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
      WRITE(IOUT,6) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
C-LFK
  160 CONTINUE
cgzh fix for 2D-problems
  165 IF(NROW.EQ.1) GO TO 175
      IF(J.NE.1) GO TO 1161
      IF(LKARR1(I,J+1,K1).EQ.0) GO TO 175
      GO TO 1163
 1161 IF(J.NE.NROW) GO TO 1162
      IF(LKARR1(I,J-1,K1).EQ.0) GO TO 175
      GO TO 1163
 1162 IF(LKARR1(I,J+1,K1).EQ.0.AND.LKARR1(I,J-1,K1).EQ.0) GO TO 175
C
C   CELL(S) LATERALLY ADJACENT TO LAKE IN Y-DIRECTION (TYPE 2) DETECTED
C
 1163 DO 170 N=1,2
      IF(N.EQ.2) GO TO 172
      IF(J.EQ.1) GO TO 170
      IF(LKARR1(I,J-1,K1).EQ.0) GO TO 170
      J2 = J - 1
      IFACE=4
      GO TO 174
  172 IF(J.EQ.NROW) GO TO 170
      IF(LKARR1(I,J+1,K1).EQ.0) GO TO 170
      J2 = J + 1
      IFACE=3
  174 M = M + 1
      IF(M.LE.MXLKND) GO TO 176
      WRITE(IOUT,149) I,J,K1
      CALL USTOP(' ')
  176 ILAKE(1,M) = K1
      ILAKE(2,M) = J
      ILAKE(3,M) = I
      ILAKE(4,M) = LKARR1(I,J2,K1)
      ILAKE(5,M) = IFACE
      BEDLAK(M) = BDLKN1(I,J,K1)
      K4 = K1 - 1
      DO 4158 K3=1,K4
      IF(LKARR1(I,J,K3).EQ.0) GO TO 4158
      GO TO 4162
 4158 CONTINUE
      BEDLAK(M) = BDLKN1(I,J,1)
 4162 CONTINUE
c-lfk      WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
C-LFK-JAN. 2013
c-lfk        WRITE(IOUT,5) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
      WRITE(IOUT,6) (ILAKE(I1,M),I1=1,5), BEDLAK(M)
C-LFK
  170 CONTINUE
  175 CONTINUE
  180 CONTINUE
      WRITE(IOUT,195) M
  195 FORMAT(/5X,'NUMBER OF LAKE-AQUIFER CELL INTERFACES = ',I5)
      LKNODE = M
C
C   SET LAKE BOTTOM ELEVATIONS
      DO 295 LK=1,NLAKES
  295 BOTTMS(LK) = 999999
C
      DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
C  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2,
C    6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      IF(NTYP.EQ.0) THEN
        LAKE = ILAKE(4,II)
Cdep  changed if statement August 24, 2009
Cdep        IF(K.GT.1) BOTLK = BOTM(I,J,LBOTM(K-1))
Cdep        IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) BOTLK = BOTM(I,J,LBOTM(K))
         IF(K.EQ.1.OR.K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) THEN
            BOTLK = BOTM(I,J,LBOTM(K))
         ELSEIF (K.EQ.0) THEN
            BOTLK = BOTM(I,J,LBOTM(1))
         ELSE
            BOTLK = BOTM(I,J,LBOTM(K-1))
         ENDIF
        IF(BOTLK.LT.BOTTMS(LAKE)) BOTTMS(LAKE) = BOTLK
      ENDIF
  350 CONTINUE
C
C-- COMPUTE AND PRINT STAGE/VOLUME TABLES WHEN MORE THAN ONE LAYER
Cdep  revised print statement to include stage/area tables
C
      IF ( IRDTAB.EQ.0 ) THEN
!      IF(NLAY.EQ.1) GO TO 1331       !RGN 5/21/12
      DO 1330 L1=1,NLAKES
      WRITE(IOUT,1306) L1
Cdep  revised print statement to include area
 1306 FORMAT(//1X,'STAGE/VOLUME RELATION FOR LAKE',I3//6X,'STAGE',
     1        8X,'VOLUME',8X,'AREA'/)
      DO  INC=1,151
        AREATABLE(INC,L1) = 0.D0
      ENDDO
      EVOL = 0.0
      GTSDPH = 40.0
      TOPMST = BOTTMS(L1)+GTSDPH
      TBELV = BOTTMS(L1)
      DO 1340 I=1,NCOL
      DO 1340 J=1,NROW
      IF(LKARR1(I,J,1).NE.L1) GO TO 1340
Cdep Revised estimate of DTHK to be thickness of top most
C     layer 6/09/2009
      IF(BOTM(I,J,0).GT.TOPMST) TOPMST = BOTM(I,J,0)
!      DTHK = BOTM(I,J,0) - BOTM(I,J,1)   RGN this was causing problems 7/8/11
!      IF (DTHK.LE.GTSDPH) THEN
!        TOPMST = BOTM(I,J,1)+DTHK
!      ELSE
!        TOPMST = BOTM(I,J,1)+GTSDPH
!      ENDIF
 1340 CONTINUE
      TBNC = (TOPMST-BOTTMS(L1))/150.0
Cdep Revised looping for computing lake stage, volume,
Cdep   and area Apr 2009.
Cdep   WRITE(IOUT,1315) TBELV, EVOL
      DO  INC=1,151
        IF (INC.GT.1) THEN
          VOLUMETABLE(INC,L1)=VOLUMETABLE(INC-1,L1)
        ENDIF
        DO I=1,NCOL
          DO J=1,NROW
            LAKEFLG = 0
            K = 1
            MOSTBOT: DO WHILE (LAKEFLG.EQ.0)
              IF(LKARR1(I,J,K).EQ.L1) THEN
                LAKEFLG = K
              ENDIF
              IF(K.EQ.NLAY)EXIT MOSTBOT
              K = K + 1
            ENDDO MOSTBOT
            IF(LAKEFLG.GT.0) THEN
              K=LAKEFLG
              FINDBOT: DO WHILE(LKARR1(I,J,K).GT.0)
                K=K+1
                IF(K.EQ.NLAY+1) EXIT
              ENDDO FINDBOT
              BOTIJ = BOTM(I,J,LBOTM(K-1))
              IF(INC.EQ.1) THEN
                IF(TBELV+1.0E-03.GT.BOTIJ) THEN
                  AREATABLE(INC,L1)=AREATABLE(INC,L1)+DELC(J)*DELR(I)
                  DEPTHTABLE(INC,L1)=TBELV
                ENDIF
              ELSE
                IF (TBELV-BOTIJ.GT.0.0) THEN
                  AREATABLE(INC,L1)=AREATABLE(INC,L1)+DELC(J)*DELR(I)
                  DEPTHTABLE(INC,L1)=TBELV
                  IF(ABS(TBELV-BOTIJ).GT.1.0E-04) THEN
                    VOLUMETABLE(INC,L1)=VOLUMETABLE(INC,L1)+
     +                                (DELC(J)*DELR(I))*TBNC
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
Cdep PRINT TABLE OF ELEVATION, VOLUME, AND AREA
        WRITE(IOUT,1315) DEPTHTABLE(INC,L1), VOLUMETABLE(INC,L1),
     +                    AREATABLE(INC,L1)
        TBELV = TBELV + TBNC
      ENDDO
 1315 FORMAT(3(1X,1PE13.5))
      WRITE(IOUT,1326)
 1326 FORMAT(120X)
Cdep  set minimum and maximum lake stages for transient simulations
      IF(ISS.EQ.0) THEN
        SSMN(L1)=BOTTMS(L1)
        SSMX(L1)=TBELV
      ENDIF
 1330 CONTINUE
 1331 CONTINUE
      ENDIF
      IF(IUNITSFR.LE.0) THEN
         NDV=0
         NTRB=0
      ENDIF
C
C
C--  READ LINKAGE PARAMETERS FOR COALESCING LAKES
C
C    FOR EACH CONNECTED LAKE SYSTEM, READ LAKE NUMBERS OF CENTER LAKES
C    AND ADJOINING LAKES AND SILL ELEVATIONS.  ENTER CARD IMAGES
C    FOR SUBLAKE SYSTEMS EVEN IF LINKED TO MAIN LAKE SYSTEM.  SYSTEMS
C    MUST BE ORDERED HIERARCHICALLY.
C
      ICMX = 0
      NCLS=0
! Read item 7
      IF(IFREFM.EQ.0) THEN
        READ(IN,'(I5)') NSLMS
      ELSE
        READ(IN,*) NSLMS
      ENDIF
      WRITE(IOUT,680) NSLMS
  680 FORMAT(/1X,'NUMBER OF CONNECTED LAKE SYSTEMS IN SIMULATION IS ',I3
     1)
      IF(NSLMS.LE.0) GO TO 760
      DO 700 IS=1,NSLMS
! Read item 8a
      IF(IFREFM.EQ.0) THEN
        READ(IN,'(16I5)',END=750) IC,(ISUB(IS,I),I=1,IC)
      ELSE
        READ(IN,*,END=750) IC,(ISUB(IS,I),I=1,IC)
      ENDIF
      IF(IC.LE.0) GO TO 750
      IF(IC.GT.ICMX) ICMX=IC
      ICS(IS)=IC
      IC1 = IC - 1
! Read item 8b
      IF(IFREFM.EQ.0) THEN
        READ(IN,'(100F10.2)') (SILLVT(IS,I),I=1,IC1)
      ELSE
        READ(IN,*) (SILLVT(IS,I),I=1,IC1)
      ENDIF
      WRITE(IOUT,18) IS, ICS(IS), ISUB(IS,1)
   18 FORMAT(/10X,'SYSTEM',I3//2X,'NUMBER OF LAKES IN SYSTEM',I5,
     1  '  CENTER LAKE NUMBER',I5//1X,'SUBLAKE NUMBER',3X,
     2  'SILL ELEVATION'/)
      DO 715 JK=2,IC
  715 WRITE(IOUT,717) ISUB(IS,JK), SILLVT(IS,JK-1)
  717 FORMAT(8X,I2,8X,F10.2)
  700 CONTINUE
  750 CONTINUE
      NCLS=IS-1
      WRITE(IOUT,751) NCLS
  751 FORMAT(/1X,'READ DATA FOR',I5,' LAKE SYSTEMS'/)
  760 CONTINUE
C
C----- READ LAKE PRECIPITATION, EVAPORATION, RUNOFF, AND WITHDRAWAL RATES.
C      IF ITMP1 LT 0, SPECIFICATIONS FROM LAST STRESS PERIOD ARE USED.
C
  800 IF(ITMP1.GE.0) GO TO 801
      WRITE(IOUT,802)
  802 FORMAT(1H0,'REUSING RECH,ET,WITHDRAWAL RATES FROM LAST STRESS PERI
     1OD'/)
      GOTO 900
  801 IF(ISS.NE.0.AND.KKPER.GT.1) WRITE(IOUT,7)
7     FORMAT(/1X,'LAKE',7X,'PRECIP',5X,'EVAP',5X,'RUNOFF',
     2     3X,'WITHDRAW',3X,'BOTTOM',5X,'AREA',5X,'SS MIN',3X,'SS MAX'
     1/90('-'))
      IF(ISS.EQ.0.OR.KKPER.EQ.1) WRITE(IOUT,77)
   77 FORMAT(/1X,'LAKE',7X,'PRECIP',5X,'EVAP',5X,'RUNOFF',
     2     3X,'WITHDRAW',3X,'BOTTOM',5X,'AREA',5X,/70('-'))
      IF (IUNITGWT.GT.0) WRITE (IOUTS,8)
 8    FORMAT (//1X,'LAKE',4X,'SOLUTE',6X,'CPPT',6X,'CRNF',6X,'CAUG'/)
! Read item 9
      DO 300 LM=1,NLAKES
      IF(IFREFM.EQ.0) THEN
        IF(ISS.NE.0.AND.KKPER.GT.1) READ(IN,'(6F10.4)') PRCPLK(LM),
     1   EVAPLK(LM),RNF(LM),WTHDRW(LM),SSMN(LM),SSMX(LM)
        IF(ISS.EQ.0.OR.KKPER.EQ.1) READ(IN,'(6F10.4)') PRCPLK(LM),
     1   EVAPLK(LM),RNF(LM),WTHDRW(LM)
      ELSE
        IF(ISS.NE.0.AND.KKPER.GT.1) READ(IN,*) PRCPLK(LM),EVAPLK(LM),
     1   RNF(LM),WTHDRW(LM),SSMN(LM),SSMX(LM)
        IF(ISS.EQ.0.OR.KKPER.EQ.1) READ(IN,*) PRCPLK(LM),EVAPLK(LM),
     1   RNF(LM),WTHDRW(LM)
      ENDIF
      IF(ISS.NE.0.AND.KKPER.GT.1) WRITE(IOUT,9) LM,PRCPLK(LM),EVAPLK(LM)
     1 ,RNF(LM),WTHDRW(LM),BOTTMS(LM),BGAREA(LM),SSMN(LM),SSMX(LM)
9     FORMAT(1X,I3,4X,1P,3E10.3,1X,5E10.3)
      IF(ISS.EQ.0.OR.KKPER.EQ.1) WRITE(IOUT,9) LM,PRCPLK(LM),EVAPLK(LM),
     1 RNF(LM),WTHDRW(LM),BOTTMS(LM),BGAREA(LM)
      IF(IUNITGWT.LE.0) GO TO 300
! Read item 9b
      DO 850 ISOL=1,NSOL
        IF(IFREFM.EQ.0) THEN
          IF(WTHDRW(LM).LT.0.0) THEN
            READ(IN,'(3F10.4)')CPPT(LM,ISOL),CRNF(LM,ISOL),CAUG(LM,ISOL)
          ELSE
            READ(IN,'(2F10.4)')CPPT(LM,ISOL),CRNF(LM,ISOL)
          ENDIF
        ELSE
          IF(WTHDRW(LM).LT.0.0) THEN
            READ(IN,*) CPPT(LM,ISOL),CRNF(LM,ISOL),CAUG(LM,ISOL)
          ELSE
            READ(IN,*) CPPT(LM,ISOL),CRNF(LM,ISOL)
          ENDIF
        ENDIF
        IF(WTHDRW(LM).LT.0.0)WRITE(IOUTS,840) LM,ISOL,
     +       CPPT(LM,ISOL),CRNF(LM,ISOL),CAUG(LM,ISOL)
        IF(WTHDRW(LM).GE.0.0)
     1  WRITE(IOUTS,841) LM,ISOL,CPPT(LM,ISOL),CRNF(LM,ISOL)
  840   FORMAT(1X,I3,6X,I3,4X,1P,3E10.2)
  841 FORMAT(1X,I3,6X,I3,4X,1P,2E10.2)
  850 CONTINUE
C      WRITE (IOUTS,'(/)')
  300 CONTINUE
      WRITE (IOUT,'(/)')
C
C------Define Initial Lake Volume & Initialize Cumulative Budget Terms
      IF(KKPER.EQ.1) THEN
!dep revised calculation of initial lake volume July 2009
        STGINIT=0.0D0
        DO 8400 LK=1,NLAKES
!dep 8400    VOL(LK)=0.0
             STGINIT=STAGES(LK)
             VOL(LK)=VOLTERP(STGINIT,LK)
             VOLINIT(LK)=VOL(LK)
 8400   CONTINUE
        DO 8450 LK=1,NLAKES
             CUMPPT(LK)=0.0
             CUMEVP(LK)=0.0
             CUMRNF(LK)=0.0
             CUMGWI(LK)=0.0
             CUMGWO(LK)=0.0
             CUMSWI(LK)=0.0
             CUMSWO(LK)=0.0
             CUMWDR(LK)=0.0
             CUMFLX(LK)=0.0
 8450   CONTINUE
            DO 8900 L=1,LKNODE
               IL=ILAKE(1,L)
               IR=ILAKE(2,L)
               IC=ILAKE(3,L)
               LAKE=ILAKE(4,L)
C------Convert ILAKE(5,L):  1 and 2 are type 1,  3 and 4 are type 2,
C        6 is type 0
               ITYPE = (ILAKE(5,L)+1)/2
               IF(ITYPE.EQ.3) ITYPE=0
               IF(ITYPE.NE.0) GO TO 8900
               IF(IL.GT.1) BOTLK = BOTM(IC,IR,LBOTM(IL-1))
               IF(IL.EQ.NLAY.AND.LKARR1(IC,IR,IL).GT.0)
     1            BOTLK = BOTM(IC,IR,LBOTM(IL))
 8900       CONTINUE
      ENDIF

 900  IF (IUNITBCF.GT.0) THEN  ! rsr, moved if block from main
        CALL SGWF2LAK7BCF7RPS()
      ELSEIF (IUNITLPF.GT.0) THEN
        CALL SGWF2LAK7LPF7RPS()
!      ELSEIF (IUNITHUF.GT.0) THEN
*        CALL SGWF2LAK7HUF7RPS()
      ELSE IF (IUNITUPW.GT.0) THEN
        CALL SGWF2LAK7UPW1RPS()
      ELSE
!        WRITE (IOUT, *) 'LAK Package requires BCF, LPF, or UPW'
        call store_error('LAK Package requires BCF, LPF, or UPW')
        CALL USTOP()
      ENDIF
      IF (IUNITSFR.GT.0) CALL SGWF2LAK7SFR7RPS()
C
C7------RETURN
      RETURN
      END SUBROUTINE GWF2LAK7RP
C
!      SUBROUTINE GWF2LAK7AD(KKPER,KKSTP,IUNITGWT,IGRID)
!C
!C------VERSION 7.1 JUNE 2006 GWF2LAK7AD; REVISED FEBRUARY 6, 2012
!C
!C     ******************************************************************
!C     ADVANCE TO NEXT TIME STEP FOR TRANSIENT LAKE SIMULATION, AND COPY
!C             INITIAL LAKE STAGES TO STGOLD FOR STEADY STATE.
!C     ******************************************************************
!C     SPECIFICATIONS:
!C     ------------------------------------------------------------------
!      USE GWFLAKMODULE, ONLY: NLAKES, LKNODE, FLOB, STAGES,
!     +                        STGNEW, STGOLD, VOLOLDD, VOLOLD, VOLINIT,
!     +                        BOTTMS, IDIV, STGOLD2, NDV
!      USE GWFSFRMODULE, ONLY: DLKSTAGE
!      USE GLOBAL,       ONLY: IOUT
!      use SimPHMFModule, only: ustop
!C     ------------------------------------------------------------------
!C     FUNCTIONS
!C     ------------------------------------------------------------------
!!*      DOUBLE PRECISION VOLTERP
!!*      EXTERNAL VOLTERP
!C     ------------------------------------------------------------------
!C
!C------SET POINTERS FOR THE CURRENT GRID.
!      CALL SGWF2LAK7PNT(IGRID)
!C
!C1 --- COPY INITIAL LAKE STAGES TO STGOLD.
!! RGN COMBINED IF AND ADDED VOLOLDD 4/17/09
!Cdep  initialized VOLINIT and VOLOLD to VOLOLDD 6/4/2009
!      DO I=1,NLAKES
!        IF(KKPER.EQ.1.AND.KKSTP.EQ.1) THEN
!          STGOLD(I)=STAGES(I)
!          VOLOLDD(I)=VOLTERP(STGOLD(I),I)
!          VOLOLD(I) = VOLOLDD(I)
!          VOLINIT(I) = VOLOLDD(I)
!          STGNEW(I)=STAGES(I)
!        ELSE
!          STGOLD2(I)=STGNEW(I)
!          STGOLD(I)=STGNEW(I)
!          VOLOLDD(I)=VOLTERP(STGOLD(I),I)
!          VOLOLD(I)=VOLOLDD(I)
!        ENDIF
!! Moved this code from 7FM  10/19/10
!         DO IDV=1,NDV
!           INODE=IDIV(I,IDV)
!           IF (INODE.GT.0) THEN
!             IF( DLKSTAGE(1,INODE).LT.DBLE(BOTTMS(I))) THEN
!               WRITE(IOUT,971)I,BOTTMS(I),
!     +                           DLKSTAGE(1,INODE),INODE
!               CALL USTOP(' ')
!             ENDIF
!           ENDIF
!         ENDDO
! ! To hear.
!      ENDDO
! 971           FORMAT(' BOTTOM ELEVATION OF LAKE ',I5,' IS ', F10.2,
!     +                 ' AND IS ABOVE OUTLET ELEVATION OF ', F10.2,
!     +                 ' FOR STREAM SEGMENT ',I5,/1X,
!     +                 ' THIS WILL CAUSE PROBLEMS IN COMPUTING LAKE',
!     +                 ' STAGE USING THE NEWTON METHOD. '/1X,
!     +                 ' ELEVATION OF STREAM OUTLET MUST BE GREATER'
!     +                 ' THAN OR EQUAL TO THE LOWEST ELEVATION OF THE',
!     +                 ' LAKE.',/1X,'*****PROGRAM STOPPING'/)
!C2 ----- IF NOT FIRST TIME STEP, OR FIRST STRESS PERIOD, UPDATE
!C           STGOLD BY STGNEW.
!! RGN MOVED TO ABOVE. STGOLD SHOULD BE UPDATED EVERY TIME STEP! 4/17/09
!!      IF (KKPER.NE.1.OR.KKSTP.NE.1) THEN
!!            DO 30 K=1,NLAKES
!!               STGOLD(K)=STGNEW(K)
!!               VOLOLD(K)=VOLTERP(STGOLD(K),K))
!!30             STGOLD2(K)=STGNEW(K)
!!      ENDIF
!C
!C-----Initialize FLOB array (stores cell by cell flux between lake and
!C                            aquifer)
!      IF (IUNITGWT.GT.0) THEN
!        DO 50 LK=1,LKNODE
! 50        FLOB(LK)=0.0
!      ENDIF
!C
!C3------RETURN
!      RETURN
!      END SUBROUTINE GWF2LAK7AD
!C
!      SUBROUTINE GWF2LAK7ST(NFLG,IGRID)
!C   ********************************************************************
!C   SET IBOUND VALUES SO THAT RECHARGE AND EVAPOTRANSPIRATION (ET) WILL
!C   BE ASSIGNED CORRECTLY UNDERNEATH DRYING LAKES (NFLG = 0), OR RESET
!C   IBOUND AFTER RECHARGE AND ET ARE COMPUTED (NFLG = 1).
!C   ********************************************************************
!C
!C   SPECIFICATIONS:
!C
!C-----------------------------------------------------------------------
!      USE GWFLAKMODULE, ONLY: LKNODE, ILAKE, STGOLD
!      USE GLOBAL,       ONLY: IBOUND, LBOTM, BOTM
!C-----------------------------------------------------------------------
!C
!C------SET POINTERS FOR THE CURRENT GRID.
!      CALL SGWF2LAK7PNT(IGRID)
!
!      IF(LKNODE.EQ.0) RETURN
!      DO 10 L=1,LKNODE
!C  Convert ILAKE(5,L):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
!      ITYPE = (ILAKE(5,L)+1)/2
!      IF(ITYPE.EQ.3) ITYPE=0
!C
!C-------ONLY CHANGE IBOUND FOR VERTICALLY ADJACENT NODE FACES
!      IF(ITYPE.NE.0) GO TO 10
!      IL = ILAKE(1,L)
!      IR = ILAKE(2,L)
!      IC = ILAKE(3,L)
!C
!C-------RESET AFTER EXECUTING RECHARGE OR ET ROUTINES
!      IF(NFLG.EQ.1) GO TO 8
!C
!C-------RESET BEFORE EXECUTING RECHARGE OR ET ROUTINES
!      IF ( IL.GT.1 ) THEN       !RGN 5/21/12 added IF test
!        IBOUND(IC,IR,IL-1) = -7
!      ELSE                      !RGN
!        IBOUND(IC,IR,IL) = -7   !RGN
!      ENDIF                    !RGN
!C
!C-------THIS IS THE CORRECT ASSIGNMENT IF PORTION OF LAKE IN COLUMN
!C       IS WET.
!      LAKE = ILAKE(4,L)
!      IF(STGOLD(LAKE).GT.BOTM(IC,IR,LBOTM(IL)-1)) GO TO 10
!C
!C-------IF PORTION OF LAKE IN NODE IS DRY, LET RECHARGE AND ET BE
!C       APPLIED TO THE AQUIFER NODE UNDERNEATH THE LAKE BY SETTING
!C       IBOUND EQUAL TO 0.
!    8 IF ( il.GT.1 ) THEN        !RGN 5/21/12 added IF test
!!     8  IBOUND(IC,IR,IL-1) = 0  !RGN
!        IBOUND(IC,IR,IL-1) = 0   !RGN
!      ELSE                       !RGN
!        IBOUND(IC,IR,IL) = 0     !RGN
!      ENDIF                     !RGN
!   10 CONTINUE
!C
!C3------RETURN
!      RETURN
!      END SUBROUTINE GWF2LAK7ST
C
      SUBROUTINE SGWF2LAK7SFR7RPS()
C
C    *******************************************************************
C--  IF STREAMS EXIST, DEFINE CONNECTIONS BETWEEN LAKES AND STREAMS
C    *******************************************************************
C
C    -------------------------------------------------------------------
C        SPECIFICATIONS:
C    -------------------------------------------------------------------
      USE GWFLAKMODULE, ONLY: NLAKES, NTRB, NDV, ITRB, IDIV, IRK
      USE GLOBAL,       ONLY: IOUT, NODES
      USE GWFSFRMODULE, ONLY: NSS, IDIVAR, IOTSG, SEG,  ISEG
      use SimPHMFModule, only: ustop
C
C-- DOUBLE CHECK SIZE OF IRK (STORED IN BUFF) vs. NLAKES
C
      IF ((NLAKES*2).GT.NODES) THEN
         WRITE (IOUT,*) '***NLAKES too large for BUFF in Subroutine GWF2
     1LAK7SFR7RPS***  STOP EXECUTION'
         CALL USTOP(' ')
      ENDIF
C
C-- INITIALIZE ARRAYS
C
      DO 55 LK=1,NLAKES
      IRK(1,LK) = 0
   55 IRK(2,LK) = 0
      NTRB = 0
      NDV = 0
C
C-- Build arrays to define lake tributary & diversion links ...
C        based on stream package input data
C
C---  Stream Inflow to Lakes
      DO 100 LSEG=1,NSS
      IF(IOTSG(LSEG).LT.0) THEN
        LAKE = -IOTSG(LSEG)
        IRK(1,LAKE) = IRK(1,LAKE) + 1
        K1 = IRK(1,LAKE)
        ITRB(LAKE,K1) = LSEG
        IF(IRK(1,LAKE).GT.NTRB) NTRB = IRK(1,LAKE)
      ENDIF
C
C---  Stream Outflow from Lakes
      IF(IDIVAR(1,LSEG).LT.0) THEN
        LAKE = -IDIVAR(1,LSEG)
        IRK(2,LAKE) = IRK(2,LAKE) + 1
        K1 = IRK(2,LAKE)
        IDIV(LAKE,K1) = LSEG
        IF(IRK(2,LAKE).GT.NDV) NDV = IRK(2,LAKE)
      ENDIF
  100 CONTINUE
C
C--  PRINT LAKE INFLOW STREAM SEGMENTS.
      WRITE(IOUT,10)
10    FORMAT(6X,'LAKE ',4X,'INFLOWING STREAM SEGMENT')
      DO 520 IK=1,NLAKES
      DO 519 JK=1,NSS
      IF(ITRB(IK,JK).LE.0) GO TO 521
  519 CONTINUE
  521 JK1 = JK - 1
      IF(JK1.GT.0) WRITE(IOUT,15) IK,(ITRB(IK,JK),JK=1,JK1)
15    FORMAT(5X,I5,14X,100I5)
  520 CONTINUE
      WRITE(IOUT,103) NTRB
103    FORMAT(/1X,'MAXIMUM NUMBER OF STREAMS INFLOWING TO A',
     1    ' LAKE IS',I5/)
C
C--  PRINT LAKE STREAM OUTFLOW SEGMENT (FROM A LAKE) NUMBERS.
C
      WRITE(IOUT,13)
13    FORMAT(6X,'LAKE ',4X,'OUTFLOWING STREAM',' SEGMENT')
      DO 600 IK=1,NLAKES
      DO 523 JK=1,NSS
      IF(IDIV(IK,JK).LE.0) GO TO 527
  523 CONTINUE
  527 JK1 = JK - 1
      IF(JK1.GT.0) WRITE(IOUT,15) IK,(IDIV(IK,JK),JK=1,JK1)
  600 CONTINUE
C
Cdep-- PRINT WARNING IF OUTFLOWING STREAM IS ASSIGNED ICALC =0.
Cdep    ADDED OCTOBER 15, 2004; DAVID PRUDIC
      DO ls = 1, NSS
        IF (IDIVAR(1,ls).LT.0) THEN
          lk = -IDIVAR(1,ls)
          IF (ISEG(1,ls).LE.0 .AND. SEG(2,ls).LE.0.0) THEN
            WRITE (IOUT, 9007) ls, lk, ISEG(1,ls), SEG(2,ls)
          ENDIF
        ENDIF
      ENDDO
      WRITE(IOUT,133) NDV
133   FORMAT(/1X,'MAXIMUM NUMBER OF STREAMS OUTFLOWING',
     1    ' FROM A LAKE IS',I5/)
 9007 FORMAT(/, ' WARNING****  OUTFLOWING STREAM SEGMENT', I6,
     +       ' FROM LAKE', I6, ' HAS AN ICALC VALUE OF', I6,
     +       ' AND FLOW INTO THE SEGMENT IS', E12.4, /,
     +       ' NO OUTFLOW FROM THE LAKE INTO ',
     +       'SEGMENT WILL BE SIMULATED', /,
     +       ' SUGGEST CHANGING ICALC TO ANOTHER OPTION')
C
C-- RETURN
      RETURN
      END SUBROUTINE SGWF2LAK7SFR7RPS
      
      SUBROUTINE SGWF2LAK7BCF7RPS()
C
C     ******************************************************************
C     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
C     THICKNESS FOR LAKES WHEN BCF PACKAGE IS USED
C     ******************************************************************
C
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
      USE GLOBAL,       ONLY: NLAY, IOUT, DELR, DELC, LAYHDT,NCOL,NROW
      USE GWFBCFMODULE, ONLY: IWDFLG, HY, CVWD, TRPY
C
      WRITE(IOUT,108)
  108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',
     1  'AQUIFER CELLS'/
     2  3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',
     3  'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/
     4  3X,'A',5X,'L',2X,'L',2X,'T',
     5  4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',
     6  'THICKNESS (L/T).)'/
     7  3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/
     8  3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',
     9  24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,
     1  'N',2X,'E',
     2  2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,
     3  'AQUIFER',2X,'COMBINED'/1X,79('_'))
C
      IWRN = 0
      IWRN1 = 0
      DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CNDFCT(II) = 0.0
C  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP = NTYP + 1
      IF(NTYP.EQ.1) THEN
C
C  Vertical Conductance
C    for vertical interface, "K" is layer below bottom of lake
C
        CNDFC1=0.0
        IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
        IF(BEDLAK(II).LE.0.0) GO TO 315
          IWRN1 = 1
        CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
        IF (IWDFLG.EQ.0) THEN
          CNDFCT(II) = CNDFC1
        ELSE
          IF(CVWD(I,J,K-1).LE.0.0.OR.CNDFC1.LE.0.0) GO TO 315
          CNDFCT(II) = 1.0/(0.5/CVWD(I,J,K-1)+1.0/CNDFC1)
        ENDIF
  315   IF (IWDFLG.EQ.0) THEN
          WRITE(IOUT,7324) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1        BEDLAK(II),CNDFC1,CNDFCT(II)
c-lfk
 7324   FORMAT(1X,5I3,2X,1P,4E10.2,10X,E11.3)
C 7324     FORMAT(1X,5I3,2X,1P,4E10.2,10X,E10.2)
        ELSE
          IF (K.GT.1) THEN
            CVWD2= 2.0*CVWD(I,J,K-1)
            WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1          BEDLAK(II),CNDFC1,CVWD2,CNDFCT(II)
          ENDIF
c-lfk
 7325   FORMAT(1X,5I3,2X,1P,5E10.2,E11.3)
c 7325     FORMAT(1X,5I3,2X,1P,6E10.2)
        ENDIF
      ELSE
C
C  Horizontal conductance
C
C  HY not read in, thus unavailable.
C
Cdep  348   IF(LAYHDT(K).EQ.0) THEN
        IF(LAYHDT(K).EQ.0) THEN
          IF(NTYP.EQ.2) CNDFCT(II) = BEDLAK(II)*DELC(J)
          IF(NTYP.EQ.3) CNDFCT(II) = BEDLAK(II)*DELR(I)
          WRITE(IOUT,7324) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1        BEDLAK(II),CNDFCT(II),CNDFCT(II)
          IWRN = 1
        ELSE
C
C  HY read in, thus available.
C
        TT = HY(I,J,K)
        IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
        IF(NTYP.EQ.3) CNDFC2 = 2.0*TRPY(K)*TT*DELR(I)/DELC(J)
        IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
        IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
        IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0)
     *         CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
        WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1    BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
        ENDIF
      ENDIF
  350 CONTINUE
C
C  WRITE WARNINGS ON LAKE/AQUIFER CONDUCTANCES, IF NECESSARY
          IF(IWRN.EQ.1.OR.IWRN1.EQ.1) WRITE(IOUT,345)
  345     FORMAT(//5X,'NOTE: INFORMATION ABOUT CALCULATED LAKE/AQUIFER C
     1ONDUCTANCES WHEN USING BCF PACKAGE FOLLOWS: '/)
          IF(IWRN.EQ.1) WRITE(IOUT,346)
  346     FORMAT(1X,'NODE(S) ADJACENT TO LAKE IN CONFINED LAYER:'/
     1    1X,'LAKE/AQUIFER CONDUCTANCES BASED SOLELY ON LAKEBED SPECIFIC
     2ATION'/)
          IF(IWRN1.EQ.1) WRITE(IOUT,347)
  347     FORMAT(1X,'IF WETDRY FLAG NOT TURNED ON, VERTICAL LEAKANCES AR
     1E NOT SAVED:'/1X,'THEREFORE, LAKE/AQUIFER CONDUCTANCES ARE BASED S
     2OLELY ON LAKEBED SPECIFICATION'/)
          IF(IWRN.EQ.1.OR.IWRN1.EQ.1) WRITE(IOUT,'(//)')
C
      RETURN
      END SUBROUTINE SGWF2LAK7BCF7RPS
C
      SUBROUTINE SGWF2LAK7LPF7RPS()
C
C     ******************************************************************
C     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
C     THICKNESS FOR LAKES WHEN LPF PACKAGE IS USED
C     ******************************************************************
C
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
      USE GLOBAL,       ONLY: NLAY, IOUT, LBOTM, LAYCBD, DELR, DELC,
     +                        BOTM
      USE GWFLPFMODULE, ONLY: CHANI, LAYVKA, VKA, VKCB, HANI, HK
C
      WRITE(IOUT,108)
  108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',
     1  'AQUIFER CELLS'/
     2  3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',
     3  'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/
     4  3X,'A',5X,'L',2X,'L',2X,'T',
     5  4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',
     6  'THICKNESS (L/T).)'/
     7  3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/
     8  3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',
     9  24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,
     1  'N',2X,'E',
     2  2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,
     3  'AQUIFER',2X,'COMBINED'/1X,79('_'))
C
      DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CAQ = 0.0
      CNDFCT(II) = 0.0
C  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP=NTYP + 1
      IF(NTYP.EQ.1) THEN
C
C  Vertical Conductance
C    for vertical interface, "K" is layer below bottom of lake
        CNDFC1=0.0
        IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
        IF(BEDLAK(II).LE.0.0) GO TO 315
        CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
        IF(LAYVKA(K).EQ.0) THEN
           VK=VKA(I,J,K)
        ELSE
           VK=HK(I,J,K)/VKA(I,J,K)
        ENDIF
c   skip if zero vk
        IF(VK.LE.0.0) GO TO 350
        BBOT=BOTM(I,J,LBOTM(K))
        TTOP=BOTM(I,J,LBOTM(K)-1)
        CAQ=VK*DELR(I)*DELC(J)/((TTOP-BBOT)*0.5)
        IF(LAYCBD(K-1).GT.0) THEN
c   skip if zero vkcb
          IF(VKCB(I,J,LAYCBD(K-1)).LE.0.0) GO TO 350
          BBOT=BOTM(I,J,LBOTM(K)-1)
          TTOP=BOTM(I,J,LBOTM(K-1))
          CCB=VKCB(I,J,LAYCBD(K-1))*DELR(I)*DELC(J)/(TTOP-BBOT)
          !include VKCB
          CAQ = 1.0/(1.0/CAQ + 1.0/CCB)
        ENDIF
        CNDFCT(II) = 1.0/(1.0/CAQ+1.0/CNDFC1)
  315   WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1             BEDLAK(II),CNDFC1,CAQ,CNDFCT(II)
      ELSE
C
C  Horizontal conductance
C
        TT = HK(I,J,K)
C X-DIRECTION
        IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
C Y-DIRECTION
        IF(NTYP.EQ.3) THEN
          IF(CHANI(K).LE.0) THEN
            KH=-CHANI(K)
            CNDFC2 = 2.0*HANI(I,J,KH)*TT*DELR(I)/DELC(J)
          ELSE
            CNDFC2 = 2.0*CHANI(K)*TT*DELR(I)/DELC(J)
          ENDIF
        ENDIF
        IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
        IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
        IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0)
     *         CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
        WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1    BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
c-lfk
 7325   FORMAT(1X,5I3,2X,1P,5E10.2,E11.3)
c 7325   FORMAT(1X,5I3,2X,1P,6E10.2)
      ENDIF
  350 CONTINUE
C
      RETURN
      END SUBROUTINE SGWF2LAK7LPF7RPS
C
      SUBROUTINE SGWF2LAK7UPW1RPS()
C
C     ******************************************************************
C     COMPUTE VERTICAL CONDUCTANCES AND HORIZONTAL CONDUCTANCES PER UNIT
C     THICKNESS FOR LAKES WHEN UPW PACKAGE IS USED
C     ******************************************************************
C
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLAKMODULE, ONLY: LKNODE, BEDLAK, LKARR1, ILAKE, CNDFCT
      USE GLOBAL,       ONLY: NLAY, IOUT, LBOTM, LAYCBD, DELR, DELC,
     +                        BOTM
      USE GWFUPWMODULE, ONLY: CHANI, LAYVKAUPW, VKAUPW, VKCB, HANI, 
     +                        HKUPW
C
      WRITE(IOUT,108)
  108 FORMAT(//9X,'C',15X,'INTERFACE CONDUCTANCES BETWEEN LAKE AND ',
     1  'AQUIFER CELLS'/
     2  3X,'L',5X,'O',10X,'(IF TYPE = 6, CONDUCTANCE (L^2/T) IS ',
     3  'BETWEEN AQUIFER CELL AND OVERLYING LAKE CELL.)',/
     4  3X,'A',5X,'L',2X,'L',2X,'T',
     5  4X,'(IF TYPE = 1 TO 4, CONDUCTANCES ARE PER UNIT SATURATED ',
     6  'THICKNESS (L/T).)'/
     7  3X,'Y',2X,'R',2X,'U',2X,'A',2X,'Y'/
     8  3X,'E',2X,'O',2X,'M',2X,'K',2X,'P',
     9  24X,'LAKEBED',6X,'C O N D U C T A N C E S'/3X,'R',2X,'W',2X,
     1  'N',2X,'E',
     2  2X,'E',5X,'DELTA Y',3X,'DELTA X',2X,'LEAKANCE',3X,'LAKEBED',3X,
     3  'AQUIFER',2X,'COMBINED'/1X,79('_'))
C
      DO 350 II=1,LKNODE
      K = ILAKE(1,II)
      J = ILAKE(2,II)
      I = ILAKE(3,II)
      CAQ = 0.0
      CNDFCT(II) = 0.0
C  Convert ILAKE(5,II):  1 and 2 are type 1,  3 and 4 are type 2, 6 is type 0
      NTYP = (ILAKE(5,II)+1)/2
      IF(NTYP.EQ.3) NTYP=0
      NTYP=NTYP + 1
      IF(NTYP.EQ.1) THEN
C
C  Vertical Conductance
C    for vertical interface, "K" is layer below bottom of lake
        CNDFC1=0.0
        IF(K.EQ.NLAY.AND.LKARR1(I,J,K).GT.0) GO TO 315
        IF(BEDLAK(II).LE.0.0) GO TO 315
        CNDFC1 = BEDLAK(II)*DELR(I)*DELC(J)
        IF(LAYVKAUPW(K).EQ.0) THEN
           VK=VKAUPW(I,J,K)
        ELSE
           VK=HKUPW(I,J,K)/VKAUPW(I,J,K)
        END IF
c   skip if zero vk
        IF(VK.LE.0.0) GO TO 350
        BBOT=BOTM(I,J,LBOTM(K))
        TTOP=BOTM(I,J,LBOTM(K)-1)
        CAQ=VK*DELR(I)*DELC(J)/((TTOP-BBOT)*0.5)
        IF(LAYCBD(K-1).GT.0) THEN
c   skip if zero vkcb
          IF(VKCB(I,J,LAYCBD(K)).LE.0.0) GO TO 350
          BBOT=BOTM(I,J,LBOTM(K)-1)
          TTOP=BOTM(I,J,LBOTM(K-1))
          CCB=VKCB(I,J,LAYCBD(K-1))*DELR(I)*DELC(J)/(TTOP-BBOT)
          !include VKCB
          CAQ = 1.0/(1.0/CAQ + 1.0/CCB)
        END IF
        CNDFCT(II) = 1.0/(1.0/CAQ+1.0/CNDFC1)
  315   WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1             BEDLAK(II),CNDFC1,CAQ,CNDFCT(II)
      ELSE
C
C  Horizontal conductance
C
        TT = HKUPW(I,J,K)
C X-DIRECTION
        IF(NTYP.EQ.2) CNDFC2 = 2.0*TT*DELC(J)/DELR(I)
C Y-DIRECTION
        IF(NTYP.EQ.3) THEN
          IF(CHANI(K).LE.0) THEN
            KH=-CHANI(K)
            CNDFC2 = 2.0*HANI(I,J,KH)*TT*DELR(I)/DELC(J)
          ELSE
            CNDFC2 = 2.0*CHANI(K)*TT*DELR(I)/DELC(J)
          END IF
        END IF
        IF(NTYP.EQ.2) CNDFC1 = BEDLAK(II)*DELC(J)
        IF(NTYP.EQ.3) CNDFC1 = BEDLAK(II)*DELR(I)
        IF (CNDFC1.GT.0.0.AND.CNDFC2.GT.0.0) 
     *         CNDFCT(II) = 1.0/(1.0/CNDFC2+1.0/CNDFC1)
        WRITE(IOUT,7325) (ILAKE(I1,II),I1=1,5),DELC(J),DELR(I),
     1    BEDLAK(II),CNDFC1,CNDFC2,CNDFCT(II)
 7325   FORMAT(1X,5I3,2X,1P,6E10.2)
      END IF
  350 CONTINUE
C
      RETURN
      END SUBROUTINE SGWF2LAK7UPW1RPS

Cdep  Added function statements to compute derivatives for Newton method
Cdep     used in solving lake stage in the FORMULATE SUBROUTINE (LAK7FM).
      DOUBLE PRECISION FUNCTION FINTERP (STAGE,LN)
Cdep&rgn  FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
C          OF LAKE STAGE TO CACULATE LAKE AREA.
C         ADDED 5/16/2006-- changed 12/2007 from "DOUBLE PRECISION FUNCTION"
C          to "FUNCTION"
      USE GWFLAKMODULE, ONLY: AREATABLE, DEPTHTABLE
      IMPLICIT NONE
      DOUBLE PRECISION STAGE, AREA, TOLF2, FOLD
      DOUBLE PRECISION a1, a2, d1, d2
      INTEGER LN, IFLG, I
      TOLF2=1.0E-7
      IF (STAGE.GT.DEPTHTABLE(151,LN))THEN
        FINTERP =  AREATABLE(151,LN)
        RETURN
      ENDIF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        a1 = AREATABLE(I,LN)
        a2 = AREATABLE(I+1,LN)
        d1 = DEPTHTABLE(I,LN)
        d2 = DEPTHTABLE(I+1,LN)
        FOLD=ABS(STAGE-d1)
        IF (FOLD .LE. TOLF2) THEN
          AREA=AREATABLE(I,LN)
          IFLG = 1
        ELSEIF (STAGE.GT.d1 .AND. STAGE.LT.d2)THEN
          AREA=((a2-a1)/(d2-d1))*STAGE+a2-((a2-a1)/(d2-d1))*d2
          IFLG = 1
        ENDIF
        I = I + 1
        IF( I.GT.150 ) THEN
          IFLG = 1
          AREA = AREATABLE(151,LN)
        ENDIF
      ENDDO
      FINTERP = AREA
      RETURN
      END FUNCTION FINTERP
      
!  RGN Added function statements to compute calculate surface area form volume
      DOUBLE PRECISION FUNCTION SURFTERP (VOLUME,LN)
C     FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
C          OF LAKE VOLUME TO CACULATE LAKE AREA.
      USE GWFLAKMODULE, ONLY: AREATABLE, VOLUMETABLE
      DOUBLE PRECISION VOLUME
      TOLF2=1.0E-7
      IF (VOLUME.GT.VOLUMETABLE(151,LN))THEN
        SURFTERP =  AREATABLE(151,LN)
        RETURN
      ENDIF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        FOLD=ABS(VOLUME-VOLUMETABLE(I,LN))
        IF (FOLD .LE. TOLF2) THEN
          AREA=AREATABLE(I,LN)
          IFLG = 1
        ELSEIF (VOLUME.GT.VOLUMETABLE(I,LN) .AND. VOLUME.LT.
     1          VOLUMETABLE(I+1,LN))THEN
          AREA=((AREATABLE(I+1,LN)-AREATABLE(I,LN))/
     1         (VOLUMETABLE(I+1,LN)- VOLUMETABLE(I,LN)))*
     2         VOLUME+AREATABLE(I+1,LN)-((AREATABLE(I+1,LN)-
     3         AREATABLE(I,LN))/(VOLUMETABLE(I+1,LN)-
     4         VOLUMETABLE(I,LN)))*VOLUMETABLE(I+1,LN)
          IFLG = 1
        ENDIF
        I = I + 1
        IF( I.GT.150 ) IFLG = 1
      ENDDO
      SURFTERP = AREA
      RETURN
      END FUNCTION SURFTERP
!
!     Interpolate lake volume as a function of lake stage
C     used in solving lake stage in the FORMULATE SUBROUTINE (LAK7FM).
      DOUBLE PRECISION FUNCTION VOLTERP (STAGE,LN)
C     FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
C          OF LAKE STAGE TO CACULATE LAKE VOLUME.
      USE GWFLAKMODULE, ONLY: VOLUMETABLE, DEPTHTABLE, AREATABLE
      IMPLICIT NONE
      INTEGER LN, IFLG, I
      DOUBLE PRECISION STAGE, VOLUME, TOLF2, FOLD
      TOLF2=1.0E-7
      IF (STAGE.GT.DEPTHTABLE(151,LN))THEN
 ! bug 5/4/09 changed FINTERP TO VOLUME
        VOLTERP =  VOLUMETABLE(151,LN)+(STAGE-DEPTHTABLE(151,LN))*
     +             AREATABLE(151,LN)
        RETURN
      ENDIF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        FOLD=ABS(STAGE-DEPTHTABLE(I,LN))
        IF (FOLD .LE. TOLF2) THEN
          VOLUME=VOLUMETABLE(I,LN)
          IFLG = 1
        ELSEIF (STAGE.GT.DEPTHTABLE(I,LN) .AND. STAGE.LT.
     1          DEPTHTABLE(I+1,LN))THEN
          VOLUME=((VOLUMETABLE(I+1,LN)-VOLUMETABLE(I,LN))/
     1         (DEPTHTABLE(I+1,LN)- DEPTHTABLE(I,LN)))*
     2         STAGE+VOLUMETABLE(I+1,LN)-((VOLUMETABLE(I+1,LN)-
     3         VOLUMETABLE(I,LN))/(DEPTHTABLE(I+1,LN)-
     4         DEPTHTABLE(I,LN)))*DEPTHTABLE(I+1,LN)
          IFLG = 1
        ENDIF
        I = I + 1
        IF( I.GT.150 ) THEN
          IFLG = 1
          VOLUME = VOLUMETABLE(151,LN)
        ENDIF
      ENDDO
      VOLTERP = VOLUME
      IF ( VOLTERP.LT.TOLF2 ) VOLTERP = TOLF2
      RETURN
      END FUNCTION VOLTERP
      
!     Interpolate lake STAGE as a function of lake VOLUME
C     used in solving lake stage in the FORMULATE SUBROUTINE (LAK7FM).
      DOUBLE PRECISION FUNCTION STGTERP (VOLUME,LN)
C     FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
C          OF LAKE VOLUME TO CACULATE LAKE STAGE.
      USE GWFLAKMODULE, ONLY: VOLUMETABLE, DEPTHTABLE,AREATABLE
      DOUBLE PRECISION VOLUME, STAGE
      TOLF2=1.0E-7
      IF (VOLUME.GT.VOLUMETABLE(151,LN))THEN
        STGTERP =  DEPTHTABLE(151,LN)+(VOLUME-VOLUMETABLE(151,LN))/
     +             AREATABLE(151,LN)
        RETURN
      ENDIF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        FOLD=ABS(VOLUME-VOLUMETABLE(I,LN))
        IF (FOLD .LE. TOLF2) THEN
          STGTERP=DEPTHTABLE(I,LN)
          IFLG = 1
        ELSEIF (VOLUME.GT.VOLUMETABLE(I,LN) .AND. VOLUME.LT.
     1          VOLUMETABLE(I+1,LN))THEN
          STGTERP=((DEPTHTABLE(I+1,LN)-DEPTHTABLE(I,LN))/
     1         (VOLUMETABLE(I+1,LN)- VOLUMETABLE(I,LN)))*
     2         VOLUME+DEPTHTABLE(I+1,LN)-((DEPTHTABLE(I+1,LN)-
     3         DEPTHTABLE(I,LN))/(VOLUMETABLE(I+1,LN)-
     4         VOLUMETABLE(I,LN)))*VOLUMETABLE(I+1,LN)
          IFLG = 1
        ENDIF
        I = I + 1
        IF( I.GT.150 ) THEN
          IFLG = 1
          STGTERP= 0.0
        ENDIF
      ENDDO
      RETURN
      END FUNCTION STGTERP
      
C------FUNCTION DERIVTERP FOR INTERPOLATING DERIVATIVE OF LAKE OUTFLOW.
      DOUBLE PRECISION FUNCTION DERIVTERP (STAGE,LSEG)
Cdep&rgn  FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
C          OF LAKE STAGE TO CACULATE LAKE OUTFLOW DERIVATIVE.
C         ADDED 5/16/2006-- changed 12/2007 from "DOUBLE PRECISION FUNCTION"
C          to "FUNCTION"
      USE GWFSFRMODULE, ONLY: DLKOTFLW, DLKSTAGE
      DOUBLE PRECISION STAGE, DEROTFLW, FOLD
      TOLF2=1.0E-7
      IF (STAGE.GT.DLKSTAGE(200,LSEG))THEN
        DERIVTERP =  DLKOTFLW(200,LSEG)
        RETURN
      ENDIF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        FOLD=ABS(STAGE-DLKSTAGE(I,LSEG))
        IF (FOLD .LE. TOLF2) THEN
          DEROTFLW=DLKOTFLW(I,LSEG)
          IFLG = 1          !rsr, changed ISFLG to IFLG
        ELSEIF (STAGE.LT.DLKSTAGE(1,LSEG)) THEN
          DEROTFLW=0.0D0
          IFLG = 1
        ELSEIF (STAGE.GT.DLKSTAGE(I,LSEG) .AND. STAGE.LT.
     1          DLKSTAGE(I+1,LSEG))THEN
          DEROTFLW=((DLKOTFLW(I+1,LSEG)-DLKOTFLW(I,LSEG))/
     1         (DLKSTAGE(I+1,LSEG)- DLKSTAGE(I,LSEG)))*
     2         STAGE+DLKOTFLW(I+1,LSEG)-((DLKOTFLW(I+1,LSEG)-
     3         DLKOTFLW(I,LSEG))/(DLKSTAGE(I+1,LSEG)-
     4         DLKSTAGE(I,LSEG)))*DLKSTAGE(I+1,LSEG)
          IFLG = 1
        ENDIF
        I = I + 1
        IF( I.GT.199) IFLG = 1
      ENDDO
      DERIVTERP = DEROTFLW
      RETURN
      END FUNCTION DERIVTERP
      
C------FUNCTION OUTFLWTERP FOR INTERPOLATING DERIVATIVE OF LAKE OUTFLOW.
      DOUBLE PRECISION FUNCTION OUTFLWTERP (STAGE,LSEG)
Cdep&rgn  FUNCTION LINEARLY INTERPOLATES BETWEEN TWO VALUES
C          OF LAKE OUTFLOW STORED IN SLKOTFLW ARRAY.
C         ADDED 5/16/2006-- changed 12/2007 from "DOUBLE PRECISION FUNCTION"
C          to "FUNCTION"
      USE GWFSFRMODULE, ONLY: SLKOTFLW, DLKSTAGE
      DOUBLE PRECISION STAGE, OUTFLOW, FOLD
      TOLF2=1.0E-9
      IF (STAGE.GT.DLKSTAGE(200,LSEG))THEN
        OUTFLWTERP =  SLKOTFLW(200,LSEG)
        RETURN
      ENDIF
      IFLG = 0
      I = 1
      DO WHILE ( IFLG.EQ.0 )
        FOLD=DABS(STAGE-DLKSTAGE(I,LSEG))
        IF (FOLD .LE. TOLF2) THEN
          OUTFLOW=SLKOTFLW(I,LSEG)
          IFLG = 1
        ELSEIF (STAGE.LT.DLKSTAGE(1,LSEG)) THEN
          OUTFLOW=0.0D0
          IFLG = 1
        ELSEIF (STAGE.GT.DLKSTAGE(I,LSEG) .AND. STAGE.LT.
     1          DLKSTAGE(I+1,LSEG))THEN
          OUTFLOW=((SLKOTFLW(I+1,LSEG)-SLKOTFLW(I,LSEG))/
     1         (DLKSTAGE(I+1,LSEG)- DLKSTAGE(I,LSEG)))*
     2         STAGE+SLKOTFLW(I+1,LSEG)-((SLKOTFLW(I+1,LSEG)-
     3         SLKOTFLW(I,LSEG))/(DLKSTAGE(I+1,LSEG)-
     4         DLKSTAGE(I,LSEG)))*DLKSTAGE(I+1,LSEG)
          IFLG = 1
        ENDIF
        I = I + 1
        IF( I.GT.199) IFLG = 1
      ENDDO
      OUTFLWTERP = OUTFLOW
      RETURN
      END FUNCTION OUTFLWTERP
C
C------FUNCTION FXLKOT_TERP FOR SMOOTHING SPECIFIED LAKE OUTFLOWS TO STREAMS.
C
      DOUBLE PRECISION FUNCTION FXLKOT_TERP(DSTAGE,Botlake,Splakout,dy)
      IMPLICIT NONE
      DOUBLE PRECISION DSTAGE,Botlake,Splakout, s, aa, ad, b, x, y, dy
      FXLKOT_TERP = 0.0D0
      s = 2.0
      x = DSTAGE-Botlake
      aa = -1.0d0/(s**2.0d0)
      ad = -2.0D0/(s**2.0d0)
      b = 2.0d0/s
      y = aa*x**2.0d0 + b*x
      dy = (ad*x + b)*Splakout
      IF ( x.LE.0.0 ) THEN
        y = 0.0D0
        dy = 0.0D0
      ELSEIF ( x-s.GT.-1.0e-14 ) THEN
        y = 1.0D0
        dy = 0.0D0
      ENDIF
      FXLKOT_TERP = y*Splakout
      END FUNCTION FXLKOT_TERP
C
      SUBROUTINE GET_FLOBOT(IC, IR, IL, ITYPE, INOFLO,CONDUC,
     1                FLOBOT,FLOBO3,FLOTOUZF,DLSTG,CLOSEZERO,H,
     2                THET1,ISS,LAKE,II,SURFDPTH,AREA,IUNITUZF,
     3                BOTLK,BOTCL,L1)
C
C     ******************************************************************
C     CALCULATE SEEPAGE BETWEEN LAKE AND GW CELLS
C     ******************************************************************
C
      USE GWFLAKMODULE
      USE GLOBAL,       ONLY: IBOUND, IOUT, LBOTM, BOTM, NLAY,LAYHDT
      USE GWFUZFMODULE, ONLY: IUZFBND,FINF,VKS
      IMPLICIT NONE
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
C     FUNCTIONS
C     -----------------------------------------------------------------
C     -----------------------------------------------------------------
C     ARGUMENTS
      DOUBLE PRECISION FLOBO3,FLOBOT,CONDUC,H,THET1,CLOSEZERO,DLSTG,
     1                 SURFDPTH,AREA,BOTLK,BOTCL,HH
      INTEGER ISS, LAKE, II, IC, IR, IL, ITYPE, IUNITUZF, L1
C     -----------------------------------------------------------------
      INTEGER ICHECK, LI, INOFLO
      DOUBLE PRECISION FLOBO1,FLOBO2,CONDMX,BOTLKUP,
     1                 BOTLKDN,FLOTOUZF,RAMPGW,RAMPSTGO,RAMPSTGN,
     2                 RAMPSTGON,HTEMP,HD,THCK,RAMPUP
C
C5C-----INITIALIZE GROUNDWATER SEEPAGE VARIABLES AND CONDUCTANCE FACTOR.
      FLOBO1 = 0.0D0
      FLOBO2 = 0.0D0
C
C6------COMPUTE SEEPAGE INTO OR OUT OF A LAKE BED NODE WHEN ITYPE=0.
C       HEAD CANNOT FALL BELOW LAKE BOTTOM
        IF (ITYPE.EQ.0) THEN
C
C6B------RAMP CONDUCTANCE ACROSS HORIZONTAL CELL FACE WHEN
C          LAKE STAGE AND GROUNDWATER HEAD NEAR LAKEBED.
          BOTLKUP = BOTLK + SURFDPTH
          BOTLKDN = BOTLK
          CONDMX = CONDUC
          HH = H
          IF ( HH.LT.BOTLKDN ) THEN
            HH = BOTLKDN
            INOFLO = 1
          ENDIF
          IF(SURFDPTH.GT.CLOSEZERO) THEN
            RAMPGW = CONDMX-(CONDMX/SURFDPTH)*
     +                           (BOTLKUP-HH)
            IF ( RAMPGW-CONDMX.GT.0.0D0 ) RAMPGW = CONDMX
            IF ( RAMPGW.LE.0.0D0 ) RAMPGW = 0.0D0
            RAMPSTGO = CONDMX-(CONDMX/SURFDPTH)*
     +                           (BOTLKUP-STGOLD(LAKE))
            IF ( RAMPSTGO-CONDMX.GT.0.0D0 ) RAMPSTGO = CONDMX
            IF ( RAMPSTGO.LE.0.0D0 ) RAMPSTGO = 0.0D0
            RAMPSTGN = CONDMX-(CONDMX/SURFDPTH)*
     +                           (BOTLKUP-STGNEW(LAKE))
            IF ( RAMPSTGN-CONDMX.GT.0.0D0 ) RAMPSTGN = CONDMX
            IF ( RAMPSTGN.LE.0.0D0 ) RAMPSTGN = 0.0D0
          ELSE
            RAMPGW=CONDMX
            RAMPSTGO=CONDMX
            RAMPSTGN=CONDMX
          ENDIF
          IF( HH-BOTLKDN.GT.CLOSEZERO ) THEN
            HTEMP = HH
          ELSE
            HTEMP=BOTLKDN
          ENDIF
C
C6C------COMPUTE LAKE SEEPAGE FOR STGOLD USING FLOBO1.
C        USE UPSTREAM WEIGHTING
          IF ( HH.LT.STGOLD(LAKE) ) THEN
            RAMPUP = RAMPSTGO
          ELSE
            RAMPUP = RAMPGW
          ENDIF
          CONDUC = RAMPUP
          IF( STGOLD(LAKE)-BOTLKDN.GT.CLOSEZERO ) THEN
            FLOBO1=CONDUC*(STGOLD(LAKE)-HTEMP)
          ELSE
            FLOBO1=CONDUC*(BOTLKDN-HTEMP)
          ENDIF
          IF ( IUNITUZF.GT.0 ) THEN
            IF ( IUZFBND(IC,IR).GT.0 )THEN
              IF (HH-BOTLK.LT.-0.5*SURFDPTH) THEN
                IF ( VKS(IC,IR)*AREA-FLOBO1.LT.CLOSEZERO )
     +                          THEN
                  FLOBO1 = VKS(IC,IR)*AREA
                ENDIF
              ENDIF
            ENDIF
          ENDIF
C
C6D------COMPUTE LAKE SEEPAGE FOR STGNEW USING FLOBO2 AND FLOBO3.
C        USE UPSTREAM WEIGHTING
          IF ( HH.LT.STGNEW(LAKE) ) THEN
            RAMPUP = RAMPSTGN
          ELSE
            RAMPUP =  RAMPGW
          ENDIF
          CONDUC = RAMPUP
          IF( STGNEW(LAKE)-BOTLKDN.GT.CLOSEZERO ) THEN
            FLOBO2 = CONDUC*(STGNEW(LAKE)-HTEMP)
            FLOBO3 = CONDUC*(STGNEW(LAKE)+DLSTG-HTEMP)
          ELSE
            FLOBO2 = CONDUC*(BOTLKDN-HTEMP)
            FLOBO3 = CONDUC*(BOTLKDN+DLSTG-HTEMP)
          ENDIF
          IF ( IUNITUZF.GT.0 ) THEN
            IF ( IUZFBND(IC,IR).GT.0 )THEN
              IF ( HH-BOTLK.LT.-0.5*SURFDPTH ) THEN
                IF ( VKS(IC,IR)*AREA-FLOBO2.LT.CLOSEZERO )
     +                           THEN
                  FLOBO2 = VKS(IC,IR)*AREA
                  FLOBO3 = VKS(IC,IR)*AREA
                ENDIF
              ENDIF
            ENDIF
          ENDIF
C
C6E------COMPUTE LAKE SEEPAGE (FLOBOT) AS A FRACTION OF FLOBO1 AND
C          FLOB02 AND FLOBO3 AS A FRACTION OF FLOBO1 AND FLOBO3.
          FLOBOT = THET1*FLOBO2 + (1.0D0-THET1)*FLOBO1
          FLOBO3 = THET1*FLOBO3 + (1.0D0-THET1)*FLOBO1
!          CONDUC = THET1*RAMPSTGN + (1.0D0-THET1)*RAMPSTGO
          IF ( IUNITUZF.GT.0 ) THEN
            IF ( IUZFBND(IC,IR).GT.0 )THEN
              IF ( HH-BOTLK.LT.-0.5*SURFDPTH ) THEN
                IF ( FLOBOT/AREA.GT.VKS(IC,IR) ) THEN
                  FLOBOT = VKS(IC,IR)*AREA
                  FLOBO3 = FLOTOUZF
                ENDIF
                FLOTOUZF = FLOBOT
                FLOBOT = 0.0D0
                CONDUC = FLOTOUZF/(STGNEW(LAKE)-BOTLK)
                FINF(IC,IR)=FLOTOUZF/AREA
              ENDIF
            ENDIF
          ENDIF
C
C7------COMPUTE SEEPAGE INTO OR OUT OF A LAKE WALL NODE
C         WHEN ITYPE=1 OR 2.
        ELSEIF ( ITYPE.EQ.1.OR.ITYPE.EQ.2 ) THEN
          IF( IBOUND(IC,IR,IL).GT.0 ) THEN
            HD = H
            IF( H.GT.BOTM(IC,IR,LBOTM(IL)-1) )
     1                         HD = BOTM(IC,IR,LBOTM(IL)-1)
C
C7B------CONDUCTANCE ACROSS VERTICAL CELL FACE DEPENDENT ON
C          SATURATED THICKNESS.
            IF ( LAYHDT(il).GT.0 ) THEN
              THCK = HD - BOTCL
            ELSE
              THCK = BOTM(IC,IR,LBOTM(IL)-1) - BOTCL
            ENDIF
            IF( THCK.LE.0.0 ) THCK = 0.0
            CONDUC = CONDUC*THCK
            IF ( H.LT.BOTM(IC,IR,LBOTM(IL)) )
     +               H = BOTM(IC,IR,LBOTM(IL))
C
C7C------COMPUTE LAKE SEEPAGE FOR STGOLD USING FLOBO1.
            IF( STGOLD(LAKE)-BOTCL.GT.CLOSEZERO ) THEN
              FLOBO1 = CONDUC*(STGOLD(LAKE)-H)
            ELSEIF ( H-BOTCL.GT.CLOSEZERO ) THEN
              FLOBO1 = CONDUC*(BOTCL-H)
            ENDIF
C
C7D------COMPUTE LAKE SEEPAGE FOR STGNEW USING FLOBO2 AND FLOBO3.
            IF( STGNEW(LAKE)-BOTCL.GT.CLOSEZERO )THEN
              FLOBO3 = CONDUC*(STGNEW(LAKE)+DLSTG-H)
              FLOBO2 = CONDUC*(STGNEW(LAKE)-H)
            ELSEIF ( H-BOTCL.GT.CLOSEZERO ) THEN
              FLOBO3 = CONDUC*(BOTCL+DLSTG-H)
              FLOBO2 = CONDUC*(BOTCL-H)
            ELSEIF ( STGNEW(LAKE)+DLSTG.GE.BOTCL )THEN
              FLOBO3 = CONDUC*(STGNEW(LAKE)+DLSTG-H)
            ENDIF
C
C7E------COMPUTE LAKE SEEPAGE (FLOBOT) AS A FRACTION OF FLOBO1 AND
C         FLOB02 AND FLOBO3 AS A FRACTION OF FLOBO1 AND FLOBO3.
            FLOBOT = THET1*FLOBO2 + (1.0D0-THET1)*FLOBO1
            FLOBO3  = THET1*FLOBO3 + (1.0D0-THET1)*FLOBO1
            SUMCNN(LAKE) = SUMCNN(LAKE) + CONDUC
          ENDIF
        ENDIF
C
C8-------SEEPAGE RATES ADDED TO MATRIX AND RESIDUAL TERMS.
C8B------COMPUTE FLWITER AND FLWITER3 DURING FIRST LOOP THROUGH
C          CALCULATIONS. NEGATIVE FLOBOT MEANS INTO LAKE
        IF ( II==1 ) THEN
          IF ( FLOBOT.LT.0.0D0 ) FLWITER(LAKE) =
     +             FLWITER(LAKE) - FLOBOT
          IF ( FLOBO3.LT.0.0D0 ) FLWITER3(LAKE) =
     +             FLWITER3(LAKE) - FLOBO3
        ENDIF
C8C------COMPUTE FLWITER AND FLOWITER3 DURING SECOND LOOP THROUGH
C          CALCULATIONS.
        IF ( II==2 ) THEN
          IF ( FLOBOT>=FLWITER(LAKE) ) THEN
            IF ( FLOBOT.GT.CLOSEZERO ) THEN
!              FLOBO2=FLWITER(LAKE)
!              FLOBOT = THET1*FLOBO2 + (1.0D0-THET1)*FLOBO1
              FLOBOT = FLWITER(LAKE)
              FLWITER(LAKE) = 0.0
              INOFLO = 1
            ENDIF
          ELSEIF ( FLOBOT.GT.CLOSEZERO )THEN
            FLWITER(LAKE) = FLWITER(LAKE) - FLOBOT
          ENDIF
          IF ( FLOTOUZF>=FLWITER(LAKE) ) THEN
            IF ( FLOTOUZF.GT.CLOSEZERO ) THEN
              FLOTOUZF=FLWITER(LAKE)
 !             FLOTOUZF = THET1*FLOTOUZF + (1.0D0-THET1)*FLOBO1
              FLWITER(LAKE) = 0.0
              INOFLO = 1
            ENDIF
          ELSEIF ( FLOTOUZF.GT.CLOSEZERO )THEN
            FLWITER(LAKE) = FLWITER(LAKE) - FLOTOUZF
          ENDIF
          IF ( FLOBO3>=FLWITER3(LAKE) ) THEN
            IF ( FLOBO3.GT.CLOSEZERO ) THEN
              FLOBO3=FLWITER3(LAKE)
 !             FLOBO3  = THET1*FLOBO3 + (1.0D0-THET1)*FLOBO1
              FLWITER3(LAKE) = 0.0
              INOFLO = 1
            ENDIF
          ELSEIF ( FLOBO3.GT.CLOSEZERO )THEN
            FLWITER3(LAKE) = FLWITER3(LAKE) - FLOBO3
          ENDIF
        ENDIF
C
C6E------COMPUTE LAKE SEEPAGE (FLOBOT) AS A FRACTION OF FLOBO1 AND
C          FLOB02 AND FLOBO3 AS A FRACTION OF FLOBO1 AND FLOBO3.
      RETURN
      END SUBROUTINE GET_FLOBOT
C
      end module GwfLakSubs
