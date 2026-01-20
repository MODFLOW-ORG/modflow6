      module GwfChdSubs

        use GWFCHDMODULE, only: SGWF2CHD7PNT, SGWF2CHD7PSV

      contains

      SUBROUTINE GWF2CHD7AR(IN,IGRID,storepar)
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR TIME-VARIANT SPECIFIED-HEAD CELLS AND
C     READ NAMED PARAMETER DEFINITIONS
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IFREFM
      USE GWFCHDMODULE,ONLY:NCHDS,MXCHD,NCHDVL,IPRCHD,NPCHD,ICHDPB,
     1                      NNPCHD,CHDAUX,CHDS
      use utl7module, only: U1DREL, U2DREL, ! UBDSV1, UBDSV2, UBDSVA,
     &                      urword, URDCOM, ! UBDSV4, UBDSVB,
     &                      ULSTRD
      CHARACTER*200 LINE
      double precision :: r
      logical, intent(in) :: storepar
C     ------------------------------------------------------------------
      ALLOCATE(NCHDS,MXCHD,NCHDVL,IPRCHD)
      ALLOCATE(NPCHD,ICHDPB,NNPCHD)
C
C1------IDENTIFY OPTION AND INITIALIZE # OF SPECIFIED-HEAD CELLS
      WRITE(IOUT,1)IN
    1 FORMAT(1X,/1X,'CHD -- TIME-VARIANT SPECIFIED-HEAD OPTION,',
     1  ' VERSION 7, 5/2/2005',/1X,'INPUT READ FROM UNIT ',I4)
      NCHDS=0
      NNPCHD=0
C
C2------READ AND PRINT MXCHD (MAXIMUM NUMBER OF SPECIFIED-HEAD
C2------CELLS TO BE SPECIFIED EACH STRESS PERIOD)
      CALL URDCOM(IN,IOUT,LINE)
      CALL UPARLSTAL(IN,IOUT,LINE,NPCHD,MXPC)
      IF(IFREFM.EQ.0) THEN
         READ(LINE,'(I10)') MXACTC
         LLOC=11
      ELSE
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXACTC,R,IOUT,IN)
      END IF
      WRITE(IOUT,3) MXACTC
    3 FORMAT(1X,'MAXIMUM OF ',I6,
     1  ' TIME-VARIANT SPECIFIED-HEAD CELLS AT ONE TIME')
C
C3------READ AUXILIARY VARIABLES AND PRINT OPTION
      ALLOCATE (CHDAUX(20))
      NAUX=0
      IPRCHD=1
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.20) THEN
            NAUX=NAUX+1
            CHDAUX(NAUX)=LINE(ISTART:ISTOP)
            WRITE(IOUT,12) CHDAUX(NAUX)
   12       FORMAT(1X,'AUXILIARY CHD VARIABLE: ',A)
         END IF
         GO TO 10
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,
     &'LISTS OF TIME-VARIANT SPECIFIED-HEAD CELLS WILL NOT BE PRINTED')
         IPRCHD = 0
         GO TO 10
      END IF
      NCHDVL=5+NAUX
C
C4------ALLOCATE SPACE FOR TIME-VARIANT SPECIFIED-HEAD LIST.
      ICHDPB=MXACTC+1
      MXCHD=MXACTC+MXPC
      ALLOCATE (CHDS(NCHDVL,MXCHD))
C
C1------READ NAMED PARAMETERS.
      WRITE(IOUT,1000) NPCHD
 1000 FORMAT(1X,//1X,I5,' TIME-VARIANT SPECIFIED-HEAD PARAMETERS')
      IF(NPCHD.GT.0) THEN
        NAUX=NCHDVL-5
        LSTSUM=ICHDPB
        DO 120 K=1,NPCHD
          LSTBEG=LSTSUM
          CALL UPARLSTRP(LSTSUM,MXCHD,IN,IOUT,IP,'CHD','CHD',1,
     &                  NUMINST,storepar)
          NLST=LSTSUM-LSTBEG
          IF (NUMINST.GT.1) NLST = NLST/NUMINST
C         ASSIGN STARTING INDEX FOR READING INSTANCES
          IF (NUMINST.EQ.0) THEN
            IB=0
          ELSE
            IB=1
          ENDIF
C         READ LIST(S) OF CELLS, PRECEDED BY INSTANCE NAME IF NUMINST>0
          LB=LSTBEG
          DO 110 I=IB,NUMINST
            IF (I.GT.0) THEN
              CALL UINSRP(I,IN,IOUT,IP,IPRCHD)
            ENDIF
            CALL ULSTRD(NLST,CHDS,LB,NCHDVL,MXCHD,0,IN,IOUT,
     &     'CHD NO.   LAYER   ROW   COL   START FACTOR      END FACTOR',
     &      CHDAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,5,IPRCHD)
            LB=LB+NLST
  110     CONTINUE
  120   CONTINUE
      END IF
C
C3------RETURN.
      CALL SGWF2CHD7PSV(IGRID)
      RETURN
      END SUBROUTINE GWF2CHD7AR

C*******************************************************************************

      SUBROUTINE GWF2CHD7RP(IN,IGRID)
C     ******************************************************************
C     READ STRESS PERIOD DATA FOR CHD
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      !use ChdModule, only: ChdType
      !use ListModule, only: ListType
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,IFREFM,IBOUND
      USE GWFCHDMODULE,ONLY:NCHDS,MXCHD,NCHDVL,IPRCHD,NPCHD,ICHDPB,
     1                      NNPCHD,CHDAUX,CHDS
      use utl7module, only: U1DREL, U2DREL, ! UBDSV1, UBDSV2, UBDSVA,
     &                      urword, URDCOM, ! UBDSV4, UBDSVB,
     &                      ULSTRD
      use SimPHMFModule, only: ustop
      ! dummy
      !type(ListType), intent(inout) :: ChdList
      ! local
      !integer :: i, j, k, nibchds
      !type(ChdType), pointer :: chd => null()
      !class(*), pointer :: obj => null()
C     ------------------------------------------------------------------
      CALL SGWF2CHD7PNT(IGRID)
      !
      ! Ned todo: If there are IBOUND-specified CHDs, include these
      !nibchds = ChdList%Count()
      !do ii=1,nibchds
      !  obj => ChdList%GetValue(ii)
      !  select type (chd => obj)
      !  type is (ChdType)
      !    j = chd%jcol
      !    i = chd%irow
      !    k = chd%klay
      !    write(iu,20)k,i,j,chd%head
      !  end select
      !enddo
C
C1------READ ITMP(FLAG TO REUSE DATA AND NUMBER OF PARAMETERS.
      IF(NPCHD.GT.0) THEN
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
C2------CALCULATE NUMBER OF AUXILIARY VALUES
      NAUX=NCHDVL-5
      IOUTU = IOUT
      IF (IPRCHD.EQ.0) IOUTU = -IOUT
C
C2------TEST ITMP
C2A-----IF ITMP<0 THEN REUSE DATA FROM LAST STRESS PERIOD
      IF(ITMP.LT.0) THEN
         WRITE(IOUT,7)
    7    FORMAT(1X,/1X,'REUSING NON-PARAMETER SPECIFIED-HEAD DATA FROM',
     1     ' LAST STRESS PERIOD')
!         nibchds = 0
      ELSE
         NNPCHD = ITMP
      END IF
C
C3------IF THERE ARE NEW NON-PARAMETER CHDS, READ THEM
      MXACTC=ICHDPB-1
      IF(ITMP.GT.0) THEN
         IF(NNPCHD.GT.MXACTC) THEN
            WRITE(IOUT,99) NNPCHD,MXACTC
   99       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE CHD CELLS (',I6,
     1                     ') IS GREATER THAN MXACTC(',I6,')')
            CALL USTOP(' ')
         END IF
         CALL ULSTRD(NNPCHD,CHDS,1,NCHDVL,MXCHD,0,IN,IOUT,
     1    'CHD NO.   LAYER   ROW   COL    START HEAD        END HEAD',
     2     CHDAUX,20,NAUX,IFREFM,NCOL,NROW,NLAY,4,5,IPRCHD)
         ! Ned todo: Add CHDs from ChdList here
      END IF
      NCHDS = NNPCHD
C
Cx------IF THERE ARE ACTIVE CHD PARAMETERS, READ THEM AND SUBSTITUTE
      CALL PRESET('CHD')
      IF(NP.GT.0) THEN
         DO 30 N=1,NP
         CALL UPARLSTSUB(IN,'CHD',IOUTU,'CHD',CHDS,NCHDVL,MXCHD,NCHDVL,
     1             MXACTC,NCHDS,4,5,
     2    'CHD NO.   LAYER   ROW   COL    START HEAD        END HEAD',
     3            CHDAUX,20,NAUX)
   30    CONTINUE
      END IF
      !
C
C4------PRINT # OF SPECIFIED-HEAD CELLS THIS STRESS PERIOD
      WRITE(IOUT,1) NCHDS
    1 FORMAT(1X,//1X,I6,' TIME-VARIANT SPECIFIED-HEAD CELLS')
C
C5------SET IBOUND NEGATIVE AT SPECIFIED-HEAD CELLS.
      DO 250 II=1,NCHDS
      IL=CHDS(1,II)
      IR=CHDS(2,II)
      IC=CHDS(3,II)
      IF(IBOUND(IC,IR,IL).GT.0) IBOUND(IC,IR,IL)=-IBOUND(IC,IR,IL)
      IF(IBOUND(IC,IR,IL).EQ.0) THEN
         WRITE(IOUT,6) IL,IR,IC
    6    FORMAT(1X,'CELL (',I3,',',I5,',',I5,') IS NO FLOW (IBOUND=0)',/
     1      1X,'NO-FLOW CELLS CANNOT BE CONVERTED TO SPECIFIED HEAD')
         CALL USTOP(' ')
      END IF
  250 CONTINUE
C
C8------RETURN
      RETURN
      END SUBROUTINE GWF2CHD7RP

      end module GwfChdSubs
