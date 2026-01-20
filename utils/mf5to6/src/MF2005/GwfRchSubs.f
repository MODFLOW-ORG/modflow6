      module GwfRchSubs
        
        use GwfRchModule, only: SGWF2RCH7PNT, SGWF2RCH7PSV
        
        private
        public :: GWF2RCH7AR, GWF2RCH7RP
        
      contains

      SUBROUTINE GWF2RCH7AR(IN,IGRID)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR RECHARGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,IFREFM
      USE GWFRCHMODULE,ONLY:NRCHOP,IRCHCB,NPRCH,IRCHPF,RECH,IRCH
      use utl7module, only: U1DREL, U2DREL, !UBDSV1, UBDSV2, UBDSVA,
     &                      urword, URDCOM, !UBDSV4, UBDSVB,
     &                      ULSTRD
      use SimPHMFModule, only: ustop
C
      CHARACTER*200 LINE
      CHARACTER*4 PTYP
      double precision :: r
C     ------------------------------------------------------------------
C
C1-------ALLOCATE SCALAR VARIABLES.
      ALLOCATE(NRCHOP,IRCHCB)
      ALLOCATE(NPRCH,IRCHPF)
C
C2------IDENTIFY PACKAGE.
      IRCHPF=0
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'RCH -- RECHARGE PACKAGE, VERSION 7, 5/2/2005',
     1' INPUT READ FROM UNIT ',I4)
C
C3------READ NRCHOP AND IRCHCB.
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARARRAL(IN,IOUT,LINE,NPRCH)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') NRCHOP,IRCHCB
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NRCHOP,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IRCHCB,R,IOUT,IN)
      END IF
C
C4------CHECK TO SEE THAT OPTION IS LEGAL.
      IF(NRCHOP.LT.1.OR.NRCHOP.GT.3) THEN
        WRITE(IOUT,8) NRCHOP
    8   FORMAT(1X,'ILLEGAL RECHARGE OPTION CODE (NRCHOP = ',I5,
     &       ') -- SIMULATION ABORTING')
        CALL USTOP(' ')
      END IF
C
C5------OPTION IS LEGAL -- PRINT OPTION CODE.
      IF(NRCHOP.EQ.1) WRITE(IOUT,201)
  201 FORMAT(1X,'OPTION 1 -- RECHARGE TO TOP LAYER')
      IF(NRCHOP.EQ.2) WRITE(IOUT,202)
  202 FORMAT(1X,'OPTION 2 -- RECHARGE TO ONE SPECIFIED NODE IN EACH',
     1     ' VERTICAL COLUMN')
      IF(NRCHOP.EQ.3) WRITE(IOUT,203)
  203 FORMAT(1X,'OPTION 3 -- RECHARGE TO HIGHEST ACTIVE NODE IN',
     1     ' EACH VERTICAL COLUMN')
C
C6------IF CELL-BY-CELL FLOWS ARE TO BE SAVED, THEN PRINT UNIT NUMBER.
      IF(IRCHCB.GT.0) WRITE(IOUT,204) IRCHCB
  204 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
C
C7------ALLOCATE SPACE FOR THE RECHARGE (RECH) AND INDICATOR (IRCH)
C7------ARRAYS.
      ALLOCATE (RECH(NCOL,NROW))
      ALLOCATE (IRCH(NCOL,NROW))
C
C8------READ NAMED PARAMETERS
      WRITE(IOUT,5) NPRCH
    5 FORMAT(1X,//1X,I5,' Recharge parameters')
      IF(NPRCH.GT.0) THEN
         DO 20 K=1,NPRCH
         CALL UPARARRRP(IN,IOUT,N,0,PTYP,1,1,0)
         IF(PTYP.NE.'RCH') THEN
            WRITE(IOUT,7)
    7       FORMAT(1X,'Parameter type must be RCH')
            CALL USTOP(' ')
         END IF
   20    CONTINUE
      END IF
C
C9------RETURN
      CALL SGWF2RCH7PSV(IGRID)
      RETURN
      END SUBROUTINE GWF2RCH7AR

!***********************************************************************

      SUBROUTINE GWF2RCH7RP(IN,IGRID)
C     ******************************************************************
C     READ RECHARGE DATA FOR STRESS PERIOD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,DELR,DELC
      USE GWFRCHMODULE,ONLY:NRCHOP,NPRCH,IRCHPF,RECH,IRCH
      use utl7module, only: U1DREL, U2DREL, !UBDSV1, UBDSV2, UBDSVA,
     &                      urword, URDCOM, !UBDSV4, UBDSVB,
     &                      ULSTRD, u2dint
      use SimPHMFModule, only: ustop
C
      CHARACTER*24 ANAME(2)
C
      DATA ANAME(1) /'    RECHARGE LAYER INDEX'/
      DATA ANAME(2) /'                RECHARGE'/
C     ------------------------------------------------------------------
C
C1------SET POINTERS FOR THE CURRENT GRID.
      CALL SGWF2RCH7PNT(IGRID)
C
C2------READ FLAGS SHOWING WHETHER DATA IS TO BE REUSED.
      IF(NRCHOP.EQ.2) THEN
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(2I10)') INRECH,INIRCH
         ELSE
            READ(IN,*) INRECH,INIRCH
         END IF
      ELSE
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10)') INRECH
         ELSE
            READ(IN,*) INRECH
         END IF
      END IF
C
C3------TEST INRECH TO SEE HOW TO DEFINE RECH.
      IF(INRECH.LT.0) THEN
C
C3A-----INRECH<0, SO REUSE RECHARGE ARRAY FROM LAST STRESS PERIOD.
        WRITE(IOUT,3)
    3   FORMAT(1X,/1X,'REUSING RECH FROM LAST STRESS PERIOD')
      ELSE
C
C3B-----INRECH=>0, SO READ RECHARGE RATE.
        IF(NPRCH.EQ.0) THEN
C
C3B1----THERE ARE NO PARAMETERS, SO READ RECH USING U2DREL.
          CALL U2DREL(RECH,ANAME(2),NROW,NCOL,0,IN,IOUT)
        ELSE
C
C3B2----DEFINE RECH USING PARAMETERS.  INRECH IS THE NUMBER OF
C3B2----PARAMETERS TO USE THIS STRESS PERIOD.
          CALL PRESET('RCH')
          WRITE(IOUT,33)
   33     FORMAT(1X,///1X,
     1      'RECH array defined by the following parameters:')
          IF(INRECH.EQ.0) THEN
              WRITE(IOUT,34)
   34         FORMAT(' ERROR: When parameters are defined for the RCH',
     &      ' Package, at least one parameter',/,' must be specified',
     &      ' each stress period -- STOP EXECUTION (GWF2RCH7RPLL)')
              CALL USTOP(' ')
            END IF
            CALL UPARARRSUB2(RECH,NCOL,NROW,0,INRECH,IN,IOUT,'RCH',
     1            ANAME(2),'RCH',IRCHPF)
        END IF
C
C4------MULTIPLY RECHARGE RATE BY CELL AREA TO GET VOLUMETRIC RATE.
        DO 50 IR=1,NROW
        DO 50 IC=1,NCOL
        RECH(IC,IR)=RECH(IC,IR)*DELR(IC)*DELC(IR)
   50   CONTINUE
      END IF
C
C5------IF NRCHOP=2 THEN A LAYER INDICATOR ARRAY IS NEEDED.  TEST INIRCH
C5------TO SEE HOW TO DEFINE IRCH.
      IF(NRCHOP.EQ.2) THEN
        IF(INIRCH.LT.0) THEN
C
C5A-----INIRCH<0, SO REUSE LAYER INDICATOR ARRAY FROM LAST STRESS PERIOD.
          WRITE(IOUT,2)
    2     FORMAT(1X,/1X,'REUSING IRCH FROM LAST STRESS PERIOD')
        ELSE
C
C5B-----INIRCH=>0, SO CALL U2DINT TO READ LAYER INDICATOR ARRAY(IRCH)
          CALL U2DINT(IRCH,ANAME(1),NROW,NCOL,0,IN,IOUT)
          DO 57 IR=1,NROW
          DO 57 IC=1,NCOL
          IF(IRCH(IC,IR).LT.1 .OR. IRCH(IC,IR).GT.NLAY) THEN
            WRITE(IOUT,56) IC,IR,IRCH(IC,IR)
   56       FORMAT(1X,/1X,'INVALID LAYER NUMBER IN IRCH FOR COLUMN',I4,
     1        '  ROW',I4,'  :',I4)
            CALL USTOP(' ')
          END IF
   57     CONTINUE
        END IF
      END IF
C
C6------RETURN
      RETURN
      END SUBROUTINE GWF2RCH7RP
      
      end module GwfRchSubs
