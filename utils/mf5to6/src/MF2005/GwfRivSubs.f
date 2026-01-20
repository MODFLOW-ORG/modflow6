      module GwfRivSubs
        
        use GwfRivModule, only: SGWF2RIV7PNT, SGWF2RIV7PSV
        
        private
        public :: GWF2RIV7AR, GWF2RIV7RP
        
      contains

      SUBROUTINE GWF2RIV7AR(IN,IGRID)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR RIVERS AND READ PARAMETER DEFINITIONS.
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
      USE GWFRIVMODULE, ONLY:NRIVER,MXRIVR,NRIVVL,IRIVCB,IPRRIV,NPRIV,
     1                       IRIVPB,NNPRIV,RIVAUX,RIVR
      use utl7module, only: U1DREL, U2DREL, !UBDSV1, UBDSV2, UBDSVA,
     &                      urword, URDCOM, !UBDSV4, UBDSVB,
     &                      ULSTRD
C
      CHARACTER*200 LINE
      double precision :: r
C     ------------------------------------------------------------------
C
C1------Allocate scalar variables, which makes it possible for multiple
C1------grids to be defined.
      ALLOCATE(NRIVER,MXRIVR,NRIVVL,IRIVCB,IPRRIV,NPRIV,IRIVPB,NNPRIV)
C
C2------IDENTIFY PACKAGE AND INITIALIZE NRIVER AND NNPRIV.
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'RIV -- RIVER PACKAGE, VERSION 7, 5/2/2005',
     1' INPUT READ FROM UNIT ',I4)
      NRIVER=0
      NNPRIV=0
C
C3------READ MAXIMUM NUMBER OF RIVER REACHES AND UNIT OR FLAG FOR
C3------CELL-BY-CELL FLOW TERMS.
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPRIV,MXPR)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(2I10)') MXACTR,IRIVCB
         LLOC=21
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTR,R,IOUT,IN)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IRIVCB,R,IOUT,IN)
      END IF
      WRITE(IOUT,3) MXACTR
    3 FORMAT(1X,'MAXIMUM OF ',I6,' ACTIVE RIVER REACHES AT ONE TIME')
      IF(IRIVCB.LT.0) WRITE(IOUT,7)
    7 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      IF(IRIVCB.GT.0) WRITE(IOUT,8) IRIVCB
    8 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
C
C4------READ AUXILIARY VARIABLES AND PRINT OPTION.
      ALLOCATE (RIVAUX(20))
      NAUX=0
      IPRRIV=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.20) THEN
            NAUX=NAUX+1
            RIVAUX(NAUX)=LINE(ISTART:ISTOP)
            WRITE(IOUT,12) RIVAUX(NAUX)
   12       FORMAT(1X,'AUXILIARY RIVER VARIABLE: ',A)
         END IF
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,'LISTS OF RIVER CELLS WILL NOT BE PRINTED')
         IPRRIV = 0
         GO TO 10
      END IF
C
C5------ALLOCATE SPACE FOR RIVER ARRAYS.
C5------FOR EACH REACH, THERE ARE SIX INPUT DATA VALUES PLUS ONE
C5------LOCATION FOR CELL-BY-CELL FLOW.
      NRIVVL=7+NAUX
      IRIVPB=MXACTR+1
      MXRIVR=MXACTR+MXPR
      ALLOCATE (RIVR(NRIVVL,MXRIVR))
C
C6------READ NAMED PARAMETERS.
      WRITE(IOUT,99) NPRIV
   99 FORMAT(1X,//1X,I5,' River parameters')
      IF(NPRIV.GT.0) THEN
        LSTSUM=IRIVPB
        DO 120 K=1,NPRIV
          LSTBEG=LSTSUM
          CALL UPARLSTRP(LSTSUM,MXRIVR,IN,IOUT,IP,'RIV','RIV',1,
     &                   NUMINST, .true.)
          NLST=LSTSUM-LSTBEG
          IF (NUMINST.EQ.0) THEN
C6A-----READ PARAMETER WITHOUT INSTANCES
            CALL ULSTRD(NLST,RIVR,LSTBEG,NRIVVL,MXRIVR,1,IN,
     &            IOUT,'REACH NO.  LAYER   ROW   COL'//
     &            '     STAGE    STRESS FACTOR     BOTTOM EL.',
     &            RIVAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRRIV)
          ELSE
C6B-----READ INSTANCES
            NINLST = NLST/NUMINST
            DO 110 I=1,NUMINST
            CALL UINSRP(I,IN,IOUT,IP,IPRRIV)
            CALL ULSTRD(NINLST,RIVR,LSTBEG,NRIVVL,MXRIVR,1,IN,
     &            IOUT,'REACH NO.  LAYER   ROW   COL'//
     &            '     STAGE    STRESS FACTOR     BOTTOM EL.',
     &            RIVAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRRIV)
            LSTBEG=LSTBEG+NINLST
  110       CONTINUE
          END IF
  120   CONTINUE
      END IF
C
C7------SAVE POINTERS TO DATA AND RETURN.
      CALL SGWF2RIV7PSV(IGRID)
      RETURN
      END SUBROUTINE GWF2RIV7AR

!***********************************************************************

      SUBROUTINE GWF2RIV7RP(IN,IGRID)
C     ******************************************************************
C     READ RIVER HEAD, CONDUCTANCE AND BOTTOM ELEVATION
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,       ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
      USE GWFRIVMODULE, ONLY:NRIVER,MXRIVR,NRIVVL,IPRRIV,NPRIV,
     1                       IRIVPB,NNPRIV,RIVAUX,RIVR
      use utl7module, only: U1DREL, U2DREL, !UBDSV1, UBDSV2, UBDSVA,
     &                      urword, URDCOM, !UBDSV4, UBDSVB,
     &                      ULSTRD
      use SimPHMFModule, only: ustop
C     ------------------------------------------------------------------
      CALL SGWF2RIV7PNT(IGRID)
C
C1------READ ITMP (NUMBER OF RIVER REACHES OR FLAG TO REUSE DATA) AND
C1------NUMBER OF PARAMETERS.
      IF(NPRIV.GT.0) THEN
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(2I10)') ITMP,NP
         ELSE
            READ(IN,*) ITMP,NP
         END IF
      ELSE
         NP=0
         IF(IFREFM.EQ.0) THEN
            READ(IN,'(I10)') ITMP
         ELSE
            READ(IN,*) ITMP
         END IF
      END IF
C
C------CALCULATE SOME CONSTANTS
      NAUX=NRIVVL-7
      IOUTU = IOUT
      IF (IPRRIV.EQ.0) IOUTU = -IOUT
C
C2------DETERMINE THE NUMBER OF NON-PARAMETER REACHES.
      IF(ITMP.LT.0) THEN
         WRITE(IOUT,7)
    7    FORMAT(1X,/1X,
     1   'REUSING NON-PARAMETER RIVER REACHES FROM LAST STRESS PERIOD')
      ELSE
         NNPRIV=ITMP
      END IF
C
C3------IF THERE ARE NEW NON-PARAMETER REACHES, READ THEM.
      MXACTR=IRIVPB-1
      IF(ITMP.GT.0) THEN
         IF(NNPRIV.GT.MXACTR) THEN
            WRITE(IOUT,99) NNPRIV,MXACTR
   99       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE REACHES (',I6,
     1                     ') IS GREATER THAN MXACTR(',I6,')')
            CALL USTOP(' ')
         END IF
         CALL ULSTRD(NNPRIV,RIVR,1,NRIVVL,MXRIVR,1,IN,IOUT,
     1          'REACH NO.  LAYER   ROW   COL'//
     2          '     STAGE      CONDUCTANCE     BOTTOM EL.',
     3          RIVAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,5,5,IPRRIV)
      END IF
      NRIVER=NNPRIV
C
C1C-----IF THERE ARE ACTIVE RIV PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('RIV')
      IF(NP.GT.0) THEN
         NREAD=NRIVVL-1
         DO 30 N=1,NP
         CALL UPARLSTSUB(IN,'RIV',IOUTU,'RIV',RIVR,NRIVVL,MXRIVR,NREAD,
     1                MXACTR,NRIVER,5,5,
     2   'REACH NO.  LAYER   ROW   COL'//
     3   '     STAGE      CONDUCTANCE     BOTTOM EL.',RIVAUX,20,NAUX)
   30    CONTINUE
      END IF
C
C3------PRINT NUMBER OF REACHES IN CURRENT STRESS PERIOD.
      WRITE (IOUT,101) NRIVER
  101 FORMAT(1X,/1X,I6,' RIVER REACHES')
C
C8------RETURN.
  260 RETURN
      END SUBROUTINE GWF2RIV7RP

      end module GwfRivSubs
