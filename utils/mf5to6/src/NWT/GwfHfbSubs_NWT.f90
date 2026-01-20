module HfBSubsNwt
  use OpenSpecModule, only: ACCESS, ACTION, FORM
  use SimPHMFModule, only: ustop
  use utl7module, only: URDCOM, URWORD, UPCASE
  private
  public :: GWF2HFB7AR, GWF2HFB7UPW, SGWF2HFB7PNT, SGWF2HFB7PSV
contains
      
      SUBROUTINE GWF2HFB7AR(INHFB,IGRID)
!C     ******************************************************************
!C     ALLOCATE ARRAY STORAGE FOR HORIZONTAL FLOW BARRIER PACKAGE
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,LAYHDT,CR,CC,BOTM,LBOTM, &
                            DELR,DELC,IOUT
      USE GWFHFBMODULE,ONLY:MXHFB,NHFB,IPRHFB,NHFBNP,NPHFB,IHFBPB,HFB
!C
      ! dummy
      INTEGER igrid, INHFB
      ! local
      INTEGER MXACTFB
      CHARACTER*16 AUX(1)
      CHARACTER*200 LINE
      double precision :: dum, r
!C     ------------------------------------------------------------------
!C
!C1------Allocate scalar data.
      ALLOCATE(MXHFB,NHFB,IPRHFB,NHFBNP,NPHFB,IHFBPB)
!C
!C2------IDENTIFY PACKAGE.
      WRITE(IOUT,1) INHFB
    1 FORMAT(1X,/1X,'HFB -- HORIZONTAL-FLOW BARRIER', &
      ' PACKAGE, NWT VERSION 1.0.9, 7/01/2014.',/, &
      '   INPUT READ FROM UNIT ',I4)
!C
!C3------READ AND PRINT NPHFB, MXFB, NHFBNP
      CALL URDCOM(INHFB,IOUT,LINE)
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPHFB,DUM,IOUT,INHFB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXFBP,DUM,IOUT,INHFB)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NHFBNP,DUM,IOUT,INHFB)
      WRITE(IOUT,500) NPHFB,MXFBP
  500 FORMAT(1X,I5,' PARAMETERS DEFINE A MAXIMUM OF ',I6, &
             ' HORIZONTAL FLOW BARRIERS')
      WRITE(IOUT,530) NHFBNP
  530 FORMAT(1X,I6,' HORIZONTAL FLOW BARRIERS NOT DEFINED BY', &
             ' PARAMETERS')
!C
!C4------LOOK FOR NOPRINT OPTION.
      IPRHFB = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INHFB)
      IF(LINE(ISTART:ISTOP).EQ.'NOPRINT') THEN
        WRITE(IOUT,3)
    3   FORMAT(1X, &
      'LISTS OF HORIZONTAL FLOW BARRIER CELLS WILL NOT BE PRINTED')
        IPRHFB = 0
      END IF
!C
!C5------CALCULATE AMOUNT OF SPACE USED BY HFB PACKAGE AND ALLOCATE HFB.
      MXACTFB = NHFBNP+MXFBP
      IHFBPB = MXACTFB + 1
      MXHFB = MXACTFB + MXFBP
      ALLOCATE (HFB(7,MXHFB))
!C
!C6------CHECK THAT THE FLOW PACKAGE IS A KIND THAT HFB CAN SUPPORT.
!C6------LAYHDT IS -1 UNLESS THE FLOW PACKAGE CHANGES IT.  IF LAYHDT
!C6------IS STILL NEGATIVE, IT IS ASSUMED THAT HFB WILL NOT WORK.
      IF (LAYHDT(1).LT.0) THEN
        WRITE(IOUT,550)
  550   FORMAT(/, &
      ' ERROR: SELECTED FLOW PACKAGE DOES NOT SUPPORT HFB PACKAGE',/, &
      ' -- STOP EXECUTION (GWF2HFB7AR)')
        CALL USTOP(' ')
      ENDIF
!C
!C7------READ PARAMETER DEFINITIONS (ITEMS 2 AND 3)
      WRITE(IOUT,600) NPHFB
  600 FORMAT(//,1X,I5,' HFB parameters')
      IF (NPHFB.GT.0) THEN
        LSTSUM = IHFBPB
        DO 20 K = 1,NPHFB
          LSTBEG = LSTSUM
          CALL UPARLSTRP(LSTSUM,MXHFB,INHFB,IOUT,IP,'HFB ','HFB ', &
                         1,NUMINST, .true.)
          IF(NUMINST.GT.0) THEN
            WRITE(IOUT,*) ' INSTANCES ARE NOT SUPPORTED FOR HFB'
            CALL USTOP(' ')
          END IF
          NLST=LSTSUM-LSTBEG
          CALL SGWF2HFB7RL(NLST,HFB,LSTBEG,MXHFB,INHFB,IOUT, &
               'BARRIER  LAYER  IROW1  ICOL1  IROW2  ICOL2     FACTOR', &
               NCOL,NROW,NLAY,IPRHFB)
          CALL SGWF2HFB7CK(LSTBEG,LSTSUM-1)
   20   CONTINUE
      ENDIF
!C
!C8------READ BARRIERS NOT DEFINED BY PARAMETERS (ITEM 4)
      NHFB = 0
      WRITE(IOUT,610) NHFBNP
  610 FORMAT(/,1X,I6,' BARRIERS NOT DEFINED BY PARAMETERS')
      IF(NHFBNP.GT.0) THEN
        CALL SGWF2HFB7RL(NHFBNP,HFB,1,MXHFB,INHFB,IOUT, &
            'BARRIER  LAYER  IROW1  ICOL1  IROW2  ICOL2    HYDCHR', &
                   NCOL,NROW,NLAY,IPRHFB)
        NHFB = NHFB + NHFBNP
        CALL SGWF2HFB7CK(1,NHFBNP)
      ENDIF
!C
!C9------SUBSTITUTE DATA FOR PARAMETERIZED BARRIERS INTO ACTIVE SECTION
!C9------OF HFB ARRAY
      IOUTU = IOUT
      IF (IPRHFB.EQ.0) IOUTU = -IOUT
      MXACTFB= IHFBPB-1
      CALL PRESET('HFB ')
      IF(NPHFB.GT.0) THEN
!C
!C10-----READ NUMBER OF ACTIVE HFB PARAMETERS (ITEM 5)
        READ(INHFB,*) NACTHFB
        IF (NACTHFB.GT.0) THEN
          DO 650 I = 1,NACTHFB
!C
!C11-----READ AND ACTIVATE AN HFB PARAMETER (ITEM 6)
            CALL SGWF2HFB7SUB(INHFB,'HFB ',IOUTU,'HFB ',HFB,7,MXHFB, &
                              MXACTFB,NHFB, &
             'BARRIER  LAYER  IROW1  ICOL1  IROW2  ICOL2     HYDCHR')
  650     CONTINUE
        ENDIF
      ENDIF
!C
!C12-----MODIFY HORIZONTAL BRANCH CONDUCTANCES FOR CONSTANT T LAYERS.
      CALL SGWF2HFB7MC()
      WRITE (IOUT,660) NHFB
  660 FORMAT(/,1X,1I6,' HFB BARRIERS')
!C
!C13-----SAVE POINTERS TO GRID AND RETURN.
      CALL SGWF2HFB7PSV(IGRID)
      RETURN
      END subroutine GWF2HFB7AR

      SUBROUTINE SGWF2HFB7MC()
!C     ******************************************************************
!C     MODIFY HORIZONTAL CONDUCTANCES (CR AND CC) FOR CONFINED LAYERS TO
!C     ACCOUNT FOR HORIZONTAL FLOW BARRIERS.  STORE UNMODIFIED HORIZONTAL
!C     CONDUCTANCES IN HFB(7,#) TO ALLOW CALCULATION OF SENSITIVITIES.
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,BOTM,LBOTM,DELR,DELC,CR,CC,LAYHDT
      USE GWFHFBMODULE,ONLY:NHFB,HFB
!C     ------------------------------------------------------------------
!C
!C1------INITIALIZE ERROR FLAG TO ZERO.
      IERFLG=0
!C
!C2----DO FOR EACH BARRIER IN RANGE.
      DO 10 II = 1,NHFB
        K = HFB(1,II)
!C
!C3------FIND ROW AND COLUMN NUMBERS OF THE TWO CELLS ON BOTH SIDES
!C3------OF THE BARRIER.
        I1 = HFB(2,II)
        J1 = HFB(3,II)
        I2 = HFB(4,II)
        J2 = HFB(5,II)
        TH0 = BOTM(J1,I1,LBOTM(K)-1) - BOTM(J1,I1,LBOTM(K))
        TH1 = BOTM(J2,I2,LBOTM(K)-1) - BOTM(J2,I2,LBOTM(K))
        THKAVG = (TH0+TH1)/2.0
        TDW = THKAVG*HFB(6,II)
!C
!C4------IF I1=I2, BARRIER IS BETWEEN TWO CELLS ON THE SAME ROW.
        IF (I1.EQ.I2) THEN
!C
!C5------IF J2-J1=1, THE TWO CELLS ARE NEXT TO ONE ANOTHER (DATA OK).
          IF ((J2-J1).EQ.1) THEN
!C
!C6------BARRIER CELLS ARE ADJACENT.
!C6------IF LAYER IS CONFINED AND BOTH CELLS ARE ACTIVE, SAVE
!C-------ORIGINAL CR FOR COMPUTING SENSITIVITIES AND MODIFY CR
            IF (LAYHDT(K).EQ.0) THEN
!C
!C7------IF CR(J1,I1,K) NOT 0, BOTH CELLS ARE ACTIVE.
              IF (CR(J1,I1,K).NE.0.) THEN
                HFB(7,II) = CR(J1,I1,K)
!C
!C8------MODIFY CR(J1,I1,K) TO ACCOUNT FOR BARRIER.
                CR(J1,I1,K) = TDW*CR(J1,I1,K)*DELC(I1)/ &
                              (TDW*DELC(I1)+CR(J1,I1,K))
              ENDIF
            ENDIF
          ENDIF
!C
!C9------IF J1=J2, BARRIER IS BETWEEN TWO CELLS ON THE SAME COLUMN.
        ELSEIF (J1.EQ.J2) THEN
!C
!C10-----IF I2-I1=1, THE TWO CELLS ARE NEXT TO ONE ANOTHER (DATA OK).
          IF ((I2-I1).EQ.1) THEN
!C
!C11-----BARRIER CELLS ARE ADJACENT.
!C11-----IF LAYER IS CONFINED AND BOTH CELLS ARE ACTIVE, SAVE
!C11-----ORIGINAL CC FOR COMPUTING SENSITIVITIES AND MODIFY CC
            IF (LAYHDT(K).EQ.0) THEN
!C
!C12-----IF CC(J1,I1,K) NOT 0, BOTH CELLS ARE ACTIVE.
              IF (CC(J1,I1,K).NE.0.) THEN
                HFB(7,II) = CC(J1,I1,K)
!C
!C13-----MODIFY CC(J1,I1,K) TO ACCOUNT FOR BARRIER.
                CC(J1,I1,K) = TDW*CC(J1,I1,K)*DELR(J1)/ &
                              (TDW*DELR(J1)+CC(J1,I1,K))
              ENDIF
            ENDIF
          ENDIF
        ENDIF
   10 CONTINUE
!C
!C14-----RETURN
      RETURN
      END subroutine SGWF2HFB7MC

      SUBROUTINE SGWF2HFB7CK(IB1,IB2)
!C     ******************************************************************
!C     CHECK HFB CELL LOCATIONS
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT
      USE GWFHFBMODULE,ONLY:HFB
!C     ------------------------------------------------------------------
!C
!C1----INITIALIZE ERROR FLAG TO ZERO.
      IERFLG=0
!C
!C2----CHECK EACH BARRIER IN RANGE.
      DO 10 II = IB1,IB2
!C
!C3------FIND ROW AND COLUMN NUMBERS OF THE TWO CELLS ON BOTH SIDES
!C3------OF THE BARRIER AND REARRANGE HFB ARRAY.
        I1 = MIN(HFB(2,II),HFB(4,II))
        J1 = MIN(HFB(3,II),HFB(5,II))
        I2 = MAX(HFB(2,II),HFB(4,II))
        J2 = MAX(HFB(3,II),HFB(5,II))
        HFB(2,II) = I1
        HFB(3,II) = J1
        HFB(4,II) = I2
        HFB(5,II) = J2
        ID = I2 - I1
        JD = J2 - J1
        IF (ID.LT.0 .OR. ID.GT.1 .OR. JD.LT.0 .OR. JD.GT.1 .OR. &
            ID.EQ.JD) THEN
!C
!C4------CELLS ARE NOT ADJACENT. PRINT ERROR MESSAGE AND SET ERROR FLAG.
   80     WRITE (IOUT,1) II-IB1+1
    1     FORMAT (1X,'ERROR DETECTED IN LOCATION DATA OF BARRIER NO. ', &
                  I6)
          IERFLG=1
        ENDIF
   10 CONTINUE
!C
!C5------HALT EXECUTION IF ERRORS ARE DETECTED.
      IF (IERFLG.EQ.1) CALL USTOP(' ')
!C
!C6------RETURN
      RETURN
      END SUBROUTINE SGWF2HFB7CK
      
      SUBROUTINE SGWF2HFB7RL(NLIST,HFB,LSTBEG,MXHFB,INPACK, &
                             IOUT,LABEL,NCOL,NROW,NLAY,IPRFLG)
!C     ******************************************************************
!C     Read and print a list of HFB barriers.
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      CHARACTER*(*) LABEL
      DIMENSION HFB(7,MXHFB)
      CHARACTER*200 LINE,FNAME
      CHARACTER*1 DASH(120)
      double precision :: factor, r, sfac
      DATA DASH/120*'-'/
      DATA NUNOPN/99/
!C     ------------------------------------------------------------------
!C
!C1------Check for and decode EXTERNAL and SFAC records.
      IN = INPACK
      ICLOSE = 0
      READ(IN,'(A)') LINE
      SFAC = 1.
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
         IN = I
         IF (IPRFLG.EQ.1) WRITE(IOUT,111) IN
  111    FORMAT(1X,'Reading list on unit ',I4)
         READ(IN,'(A)') LINE
      ELSE IF(LINE(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME = LINE(ISTART:ISTOP)
         IN = NUNOPN
         IF (IPRFLG.EQ.1) WRITE(IOUT,115) IN,FNAME
  115    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         OPEN(UNIT=IN,FILE=FNAME,ACTION=ACTION(1))
         ICLOSE = 1
         READ(IN,'(A)') LINE
      END IF
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SFAC') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SFAC,IOUT,IN)
         IF (IPRFLG.EQ.1) THEN
            WRITE(IOUT,116) SFAC
  116       FORMAT(1X,'LIST SCALING FACTOR= ',1PG12.5)
         ENDIF
         READ(IN,'(A)') LINE
      END IF
!C
!C2------Define label for printout.
      NBUF = LEN(LABEL)+3
      IF (IPRFLG.EQ.1) THEN
         WRITE(IOUT,103) LABEL
         WRITE(IOUT,104) (DASH(J),J=1,NBUF)
  103    FORMAT(1X,/1X,A)
  104    FORMAT(1X,400A)
      ENDIF
!C
!C3------Loop through the number of cells to read.
      N = NLIST+LSTBEG-1
      DO 250 II=LSTBEG,N
!C
!C4------Read a line into the buffer.  (The first line has already been read
!C4------in order to scan for EXTERNAL and SFAC records.)
      IF(II.NE.LSTBEG) READ(IN,'(A)') LINE
!C
!C5------Read the non-optional values from a line.
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,K,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I1,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J1,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I2,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J2,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,FACTOR,IOUT,IN)
      HFB(1,II) = K
      HFB(2,II) = I1
      HFB(3,II) = J1
      HFB(4,II) = I2
      HFB(5,II) = J2
      HFB(6,II) = FACTOR*SFAC
      HFB(7,II) = 0.0
!C
!C6------Write the values that were read.
      NN = II-LSTBEG+1
      IF (IPRFLG.EQ.1) WRITE(IOUT,205) NN,K,I1,J1,I2,J2,HFB(6,II)
205   FORMAT(1X,I6,2X,I5,1X,4(2X,I5),2X,1PG11.4)
!C
!C7------Check for illegal grid location.
      IF(K.LT.1 .OR. K.GT.NLAY) THEN
         WRITE(IOUT,*) ' Layer number in list is outside of the grid'
         CALL USTOP(' ')
      END IF
      IF(I1.LT.1 .OR. I1.GT.NROW .OR. I2.LT.1 .OR. I2.GT.NROW) THEN
         WRITE(IOUT,*) ' Row number in list is outside of the grid'
         CALL USTOP(' ')
      END IF
      IF(J1.LT.1 .OR. J1.GT.NCOL .OR. J2.LT.1 .OR. J2.GT.NCOL) THEN
         WRITE(IOUT,*) ' Column number in list is outside of the grid'
         CALL USTOP(' ')
      END IF
  250 CONTINUE
      IF(ICLOSE.NE.0) CLOSE(UNIT=IN)
!C
!C8------Return.
      RETURN
      END SUBROUTINE SGWF2HFB7RL
      
      SUBROUTINE SGWF2HFB7SUB(IN,PACK,IOUTU,PTYP,HFB,LSTVL,MXHFB, &
                      MXACTFB,NHFB,LABEL)
!C     ******************************************************************
!C     Read a parameter name, look it up in the list of parameters,
!C     and substitute values into active part of HFB array.
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE PARAMMODULE
      CHARACTER*(*) PACK,PTYP
      DIMENSION HFB(LSTVL,MXHFB)
      CHARACTER*(*) LABEL
      CHARACTER*200 LINE
      CHARACTER*10 CTMP1,CTMP2,CTMP3,CTMP4
      double precision :: rdum
!C     ------------------------------------------------------------------
!C
!C1------The Listing File file unit is the absolute value of IOUTU.  
!C1------Read the parameter name.
      IOUT = ABS(IOUTU)
      READ(IN,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
      WRITE(IOUT,1) LINE(ISTART:ISTOP)
    1 FORMAT(/,' Parameter:  ',A)
      IF(LINE(ISTART:ISTOP).EQ.' ') THEN
        WRITE(IOUT,*) ' Blank parameter name in the ',PACK,' file.'
        CALL USTOP(' ')
      END IF
!C
!C2------Find the parameter in the list of parameters.
      CTMP1=LINE(ISTART:ISTOP)
      CALL UPCASE(CTMP1)
      DO 100 IP=1,IPSUM
        CTMP2=PARNAM(IP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
          IF(PARTYP(IP).NE.PTYP) THEN
            WRITE(IOUT,11) PARNAM(IP),PARTYP(IP),PACK,PTYP
   11       FORMAT(1X,'Parameter type conflict:',/ &
              1X,'Named parameter:',A,' was defined as type:',A,/ &
              1X,'However, this parameter is used in the ',A, &
                ' file, so it should be type:',A)
            CALL USTOP(' ')
          END IF
!C
!C3------Set indices to point to the barriers that correspond to the
!C3------specified parameter.
          NLST=IPLOC(2,IP)-IPLOC(1,IP)+1
          NI=1
!C
!C4------Check that the parameter is not already active.
          IF (IACTIVE(IP).GT.0) THEN
            WRITE(IOUT,73) PARNAM(IP)
   73       FORMAT(/,1X,'*** ERROR: PARAMETER "',A, &
                '" HAS ALREADY BEEN ACTIVATED THIS STRESS PERIOD',/, &
                ' -- STOP EXECUTION (UPARLSTSUB)')
            CALL USTOP(' ')
          ENDIF
!C
!C5------Set the active flag.
          IACTIVE(IP)=NI
!C
!C6------Accumulate the total number of active barriers in the list.
          NHFB=NHFB+NLST
          IF(NHFB.GT.MXACTFB) THEN
            WRITE(IOUT,83) NHFB,MXACTFB
   83       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE LIST ENTRIES (',I6, &
             ')',/1X,'IS GREATER THAN THE MAXIMUM ALLOWED (',I6,')')
            CALL USTOP(' ')
          END IF
!C
!C7------Write label for barrier values if IOUTU is positive.
          IF(IOUTU.GT.0) THEN
             WRITE(IOUT,'(1X,A)') LABEL
             WRITE(IOUT,84)
   84        FORMAT(1X,56('-'))
           END IF
!C
!C8------Copy the values from the paramter location into the front part
!C8------of the list where the currently active list is kept.
          DO 90 I=1,NLST
            II=NHFB-NLST+I
            III=I-1+IPLOC(1,IP)+(NI-1)*NLST
            DO 85 J=1,7
              HFB(J,II)=HFB(J,III)
   85       CONTINUE
!C
!C8A-----Scale HYDCHR by the parameter value.
            HFB(6,II)=HFB(6,II)*B(IP)
            IL=HFB(1,II)
            IR1=HFB(2,II)
            IC1=HFB(3,II)
            IR2=HFB(4,II)
            IC2=HFB(5,II)
            IF(IOUTU.GT.0) WRITE(IOUT,89) II,IL,IR1,IC1,IR2,IC2, &
                HFB(6,II)
   89       FORMAT(1X,I6,2X,I5,1X,4(2X,I5),2X,1PG11.4)
   90     CONTINUE
!C
!C8B------After moving the data, return.
          RETURN
        END IF
  100 CONTINUE
!C
!C9------All parameter names have been checked without finding the
!C9------parameter. Write an error message and stop.
      WRITE(IOUT,*) ' The ',PACK, &
         ' file specifies an undefined parameter:',LINE(ISTART:ISTOP)
      CALL USTOP(' ')
!C
      END SUBROUTINE SGWF2HFB7SUB

  SUBROUTINE SGWF2HFB7PSV(IGRID)
!     Save pointers to HFB data for a grid.
  USE GWFHFBMODULE
!
    GWFHFBDAT(IGRID)%MXHFB=>MXHFB
    GWFHFBDAT(IGRID)%NHFB=>NHFB
    GWFHFBDAT(IGRID)%IPRHFB=>IPRHFB
    GWFHFBDAT(IGRID)%NHFBNP=>NHFBNP
    GWFHFBDAT(IGRID)%NPHFB=>NPHFB
    GWFHFBDAT(IGRID)%IHFBPB=>IHFBPB
    GWFHFBDAT(IGRID)%HFB=>HFB
!
  RETURN
  END SUBROUTINE SGWF2HFB7PSV

      SUBROUTINE GWF2HFB7UPW(IGRID)
!C     ******************************************************************
!C     MODIFY HORIZONTAL BRANCH CONDUCTANCES IN VARIABLE-TRANSMISSIVITY
!C     LAYERS TO ACCOUNT FOR HORIZONTAL FLOW BARRIERS IN MODFLOW-NWT. 
!C     STORE UNMODIFIED HORIZONTAL CONDUCTANCE IN HFB(7,#) TO ALLOW 
!C     CALCULATION OF SENSITIVITIES.
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,HNEW,LAYHDT,CR,CC,BOTM,LBOTM, &
                            DELR,DELC
      USE GWFHFBMODULE,ONLY:NHFB,HFB
!C     ------------------------------------------------------------------
!C
!C1------Set pointers to the specified grid.
      CALL SGWF2HFB7PNT(IGRID)
!C
!C2------FOR EACH BARRIER, MODIFY HORIZONTAL BRANCH CONDUCTANCES IF LAYER
!C2------IS CONVERTIBLE.
      DO 10 II=1,NHFB
        K = HFB(1,II)
!C
!C3-----IF LAYHDT=0, THICKNESS AND CONDUCTANCE DO NOT VARY, AND
!C3-----MODIFICATION OF CONDUCTANCE DUE TO BARRIER WAS DONE IN
!C3-----SGWF1HFBMC
        IF (LAYHDT(K).GT.0) THEN
!C
!C4------CELL (J1,I1,K) IS THE ONE WHOSE HORIZONTAL BRANCH
!C4------CONDUCTANCES ARE TO BE MODIFIED.
          I1 = HFB(2,II)
          J1 = HFB(3,II)
!C
!C5------CELL (J2,I2,K) IS THE CELL NEXT TO CELL (J1,I1,K) AND
!C5------SEPARATED FROM IT BY THE BARRIER.
          I2 = HFB(4,II)
          J2 = HFB(5,II)
          HCDW = HFB(6,II)
!C
!C6------IF I1=I2, MODIFY HORIZONTAL BRANCH CONDUCTANCES ALONG ROW
!C6------DIRECTION.
          IF (I1.EQ.I2) THEN
!C
!C7------IF CR(J1,I1,K) NOT = 0, CELLS ON EITHER SIDE OF BARRIER ARE
!C7------ACTIVE
            IF (CR(J1,I1,K).NE.0.) THEN
!C
!C
!C9------STORE UNMODIFIED CR FOR CALCULATING SENSITIVITIES
              HFB(7,II) = CR(J1,I1,K)
!C
!C10-----MODIFY CR(J1,I1,K) TO ACCOUNT FOR BARRIER.
              TDW = HCDW
              CR(J1,I1,K) = TDW*CR(J1,I1,K)*DELC(I1)/ &
                            (TDW*DELC(I1)+CR(J1,I1,K))
            ENDIF
!C
!C11-----CASE OF J1=J2. MODIFY HORIZONTAL BRANCH CONDUCTANCES ALONG
!C11-----COLUMN DIRECTION.
          ELSE
!C
!C12-----IF CC(J1,I1,K) NOT = 0, CELLS ON EITHER SIDE OF BARRIER ARE
!C12-----ACTIVE
            IF (CC(J1,I1,K).NE.0.) THEN
!C
!C14-----STORE UNMODIFIED CC FOR CALCULATING SENSITIVITIES.
              HFB(7,II) = CC(J1,I1,K)
!C
!C15-----MODIFY CC(J1,I1,K) TO ACCOUNT FOR BARRIER.
              TDW = HCDW
              CC(J1,I1,K) = TDW*CC(J1,I1,K)*DELR(J1)/ &
                            (TDW*DELR(J1)+CC(J1,I1,K))
            ENDIF
          ENDIF
        ENDIF
   10 CONTINUE
!C
!C16-----RETURN
      RETURN
      END SUBROUTINE GWF2HFB7UPW

  SUBROUTINE SGWF2HFB7PNT(IGRID)
!C  Set pointers to HFB data for a grid.
  USE GWFHFBMODULE
!C
    MXHFB=>GWFHFBDAT(IGRID)%MXHFB
    NHFB=>GWFHFBDAT(IGRID)%NHFB
    IPRHFB=>GWFHFBDAT(IGRID)%IPRHFB
    NHFBNP=>GWFHFBDAT(IGRID)%NHFBNP
    NPHFB=>GWFHFBDAT(IGRID)%NPHFB
    IHFBPB=>GWFHFBDAT(IGRID)%IHFBPB
    HFB=>GWFHFBDAT(IGRID)%HFB
    NumHFBs => GWFHFBDAT(igrid)%NumHFBs
!C
  RETURN
  END SUBROUTINE SGWF2HFB7PNT

end module HfBSubsNwt
