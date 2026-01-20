      SUBROUTINE UPARARRAL(IN,IOUT,LINE,NP)
C     ******************************************************************
C     Setup array parameter definition for a package.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      CHARACTER*(*) LINE
      double precision :: r
C     ------------------------------------------------------------------
C
C  If NP has not already been defined, decode PARAMETER definitions if
C  they exist
      IF(IN.GT.0) THEN
         NP=0
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(LINE(ISTART:ISTOP).EQ.'PARAMETER') THEN
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NP,R,IOUT,IN)
            READ(IN,'(A)') LINE
         END IF
      END IF
C
C  Process the parameter information
      IF(NP.GT.0) THEN
         WRITE(IOUT,31) NP
   31    FORMAT(1X,I5,' Named Parameters     ')
      ELSE
         NP=0
         WRITE(IOUT,'(A)') ' No named parameters'
      END IF
C
      RETURN
      END
      SUBROUTINE UPARARRRP(IN,IOUT,NP,ILFLG,PTYP,ITERP,ITVP,IACT)
C     ******************************************************************
C     Read and store array parameter definition information for one
C     parameter.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      CHARACTER*(*) PTYP
      CHARACTER*200 LINE
      CHARACTER*10 PN,CTMP1,CTMP2
      double precision :: r, pv
C     ------------------------------------------------------------------
C
C1------Read a parameter definition line and decode the parameter name,
C1------type, and value
      READ(IN,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
      PN=LINE(ISTART:ISTOP)
      CTMP1=PN
      CALL UPCASE(CTMP1)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      PTYP=LINE(ISTART:ISTOP)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,PV,IOUT,IN)
C
C2------Look for the parameter name in the parameter list
      DO 10 NP=1,MXPAR
        CTMP2=PARNAM(NP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
C
C2A-----If found, determine if it is an illegal duplicate or if it was
C         predefined.
          IF(PARTYP(NP).NE.' ' .AND. ITERP.EQ.1) THEN
C           Illegal duplicate
            WRITE(IOUT,110) CTMP1
  110       FORMAT(' Duplicate parameter name: ',A)
            CALL USTOP(' ')
          END IF
C         Parameter was predefined -- leave its value alone
C         (i.e. ignore PV).
          GO TO 100
        ELSE IF(PARNAM(NP).EQ.' ') THEN
C         Parameter was not found in the list, so it is a new
C         definition. Put values in the list.
          PARNAM(NP)=PN
          B(NP)=PV
          IPSUM=IPSUM+1
          GO TO 100
        END IF
10    CONTINUE
C
C2B-----Entire parameter list has been searched without finding
C2B-----a blank entry for the new parameter.  Too many parameters
      WRITE(IOUT,11)
   11 FORMAT(1X,'The number of parameters has exceeded the maximum')
      CALL USTOP(' ')
C
C3------Parameter is a new parameter or it was prefined in the
C3------Parameter Value file.  Get the number of clusters.
  100 PARTYP(NP)=PTYP
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NCLU,R,IOUT,IN)
      IF(NCLU.LE.0) THEN
        WRITE(IOUT,104) PN
  104   FORMAT(' ERROR:  DEFINITION FOR PARAMETER "',A,'"',
     &   ' INCLUDES NO CLUSTERS',/,'   -- STOP EXECUTION (UPARARRRP)')
        CALL USTOP(' ')
      ENDIF
      IF(ITERP.EQ.1) THEN
        NUMINST=0
        IF (ITVP.GT.0) THEN
C
C4------CHECK FOR MULTIPLE INSTANCES.
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
          IF (LINE(ISTART:ISTOP).EQ.'INSTANCES') THEN
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUMINST,R,IOUT,IN)
            IF (NUMINST.LT.1) THEN
              WRITE(IOUT,12) PARNAM(NP),PTYP
   12         FORMAT(/,1X,'*** ERROR: NUMINST SPECIFIED LESS THAN 1',
     &            ' FOR PARAMETER "',A,'"',/,12X,'OF TYPE "',A,
     &            '" -- STOP EXECUTION (UPARARRRP)')
              CALL USTOP(' ')
            ENDIF
          ENDIF
        ENDIF
C
C5------SET IPLOC VALUES.
        IPLOC(1,NP)=ICLSUM+1
        NI=MAX(1,NUMINST)
        ICLSUM=ICLSUM+NCLU*NI
        IPLOC(2,NP)=ICLSUM
        IPLOC(3,NP)=NUMINST
        IPLOC(4,NP)=INAMLOC
        INAMLOC=INAMLOC+NUMINST
C
C6------MAKE SURE THAT THE MAXIMUM NUMBER OF CLUSTERS IN IPCLST IS
C6------NOT EXCEEDED.
        IF(IPLOC(2,NP).GT.MXCLST) THEN
          WRITE(IOUT,117) IPLOC(2,NP),MXCLST
  117     FORMAT(1X,I5,
     1   ' CLUSTERS WERE SPECIFIED, BUT THERE IS SPACE FOR ONLY',I5)
          WRITE(IOUT,*) NP,NCLU
          WRITE(IOUT,'(A)') PARNAM(NP)
          WRITE(IOUT,'(4I10)') IPLOC
          CALL USTOP(' ')
        END IF
        WRITE(IOUT,121) PARNAM(NP),PARTYP(NP),NCLU
  121   FORMAT(1X/,1X,'PARAMETER NAME:',A,'   TYPE:',A,'   CLUSTERS:',
     &         I4)
        WRITE(IOUT,122) PV
  122   FORMAT(1X,'Parameter value from package file is: ',1PG13.5)
        IF(B(NP).NE.PV) THEN
          WRITE(IOUT,123) B(NP)
  123     FORMAT(1X,'This value has been changed to:',7X,1PG13.5,
     &        ', as read from',/,' the Parameter Value file')
        END IF
C
C7------MAKE SURE THE MAXIMUM NUMBER OF INSTANCES IS NOT EXCEEDED.
        IF(NUMINST.GT.0) THEN
          WRITE(IOUT,124)NUMINST
  124      FORMAT(3X,'NUMBER OF INSTANCES: ',I4)
          IF((INAMLOC-1).GT.MXINST) THEN
            WRITE(IOUT,125)INAMLOC-1,MXINST
  125       FORMAT(1X,'EXCEEDED THE MAXIMUM NUMBER OF INSTANCES:'/
     &       1X,I5,' instances have been specified'/
     &       1X,'The maximum number of instances is',I5)
            CALL USTOP(' ')
          ENDIF
        ENDIF
      ELSE
        NUMINST=IPLOC(3,NP)
      ENDIF
      IACTIVE(NP)=IACT
C
C8------Process clusters for each instance.
      IF(NUMINST.EQ.0) THEN
        IB=0
      ELSE
        IB=1
      ENDIF
      I=IPLOC(1,NP)-1
      DO 210 INST=IB,NUMINST
        IF(NUMINST.GT.0) CALL UINSRP(INST,IN,IOUT,NP,ITERP)
C
C9------Read and process clusters.
        DO 200 KK=1,NCLU
          I=I+1
          READ(IN,'(A)') LINE
          IF(ITERP.EQ.1) THEN
            LLOC=1
            IF(ILFLG.NE.0) THEN
C
C9A-----Get layer number for cluster
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPCLST(1,I),R,IOUT,
     &                    IN)
            ELSE
              IPCLST(1,I)=0
            END IF
C
C9B-----Get multiplier and zone array names.
            CALL URWORD(LINE,LLOC,IM1,IM2,0,N,R,IOUT,IN)
            CALL URWORD(LINE,LLOC,IZ1,IZ2,0,N,R,IOUT,IN)
C
C9C-----Get zone numbers.
            DO 30 J=5,14
              CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPCLST(J,I),R,-1,IN)
              IF(IPCLST(J,I).EQ.0) THEN
                IPCLST(4,I)=J-1
                GO TO 32
              END IF
   30       CONTINUE
            IPCLST(4,I)=14
   32       CONTINUE
            IF(ILFLG.NE.0) THEN
              WRITE(IOUT,36) IPCLST(1,I),LINE(IM1:IM2),LINE(IZ1:IZ2)
   36         FORMAT(16X,'LAYER: ',I3,'    MULTIPLIER ARRAY: ',A,
     2        '    ZONE ARRAY: ',A)
            ELSE
              WRITE(IOUT,37) LINE(IM1:IM2),LINE(IZ1:IZ2)
   37         FORMAT(16X,'MULTIPLIER ARRAY: ',A,'    ZONE ARRAY: ',A)
            END IF
C
C9D-----Find the multiplier array number.
            CTMP1=LINE(IM1:IM2)
            CALL UPCASE(CTMP1)
            IF(CTMP1.EQ.'NONE') THEN
              IPCLST(2,I)=0
            ELSE
              IF(NMLTAR.GT.0) THEN
              DO 40 J=1,NMLTAR
                CTMP2=MLTNAM(J)
                CALL UPCASE(CTMP2)
                IF(CTMP1.EQ.CTMP2) GO TO 45
   40           CONTINUE
              END IF
              WRITE(IOUT,'(A)') ' Multiplier array has not been defined'
              CALL USTOP(' ')
   45         IPCLST(2,I)=J
            END IF
C
C9E-----Find the zone array number.
            CTMP1=LINE(IZ1:IZ2)
            CALL UPCASE(CTMP1)
            IF(CTMP1.EQ.'ALL') THEN
              IPCLST(3,I)=0
            ELSE
              IF(IPCLST(4,I).EQ.4) THEN
                WRITE(IOUT,47)
   47           FORMAT(1X,
     1          'There were no zone values specified in the cluster',/
     2          1X,'At least one zone must be specified')
                CALL USTOP(' ')
              END IF
              WRITE(IOUT,48) (IPCLST(J,I),J=5,IPCLST(4,I))
   48         FORMAT(1X,'               ZONE VALUES:',10I5)
              IF(NZONAR.GT.0) THEN
                DO 50 J=1,NZONAR
                  CTMP2=ZONNAM(J)
                  CALL UPCASE(CTMP2)
                  IF(CTMP1.EQ.CTMP2) GO TO 55
   50           CONTINUE
              END IF
              WRITE(IOUT,'(A)') ' Zone array has not been defined'
              CALL USTOP(' ')
   55         IPCLST(3,I)=J
            END IF
          ENDIF
C
  200   CONTINUE
  210 CONTINUE
C
C10-----RETURN.
      RETURN
      END
      SUBROUTINE UPARARRSUB1(ZZ,NCOL,NROW,ILAY,PTYP,IOUT,ANAME,IPF)
C     ******************************************************************
C     Substitute parameter-based values into a 2-D array based on a
C     parameter type.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      DIMENSION ZZ(NCOL,NROW)
      CHARACTER*(*) PTYP
      CHARACTER*24 ANAME
      double precision :: zz
C     ------------------------------------------------------------------
C
C1------Set initialization flag to cause USUB2D to initialze ZZ to 0.
C1------Write a header above the list of parameters that define ZZ.
      INIT=1
      WRITE(IOUT,11) ANAME
   11 FORMAT(1X,/,1X,A,' is defined by the following parameters:')
C
C2------Loop through each parameter looking for the specified file type.
      DO 100 IP=1,IPSUM
C
C2A-----Stop looping if the end of the parameter list is found.
      IF(PARNAM(IP).EQ.' ') GO TO 200
C
C2B-----Check for the specified parameter type.
      IF(PARTYP(IP).EQ.PTYP) THEN
C
C2C-----Loop through each cluster definition for layers that match the
C2C-----specified layer.
         II=IP
         CALL USUB2D(ZZ,NCOL,NROW,II,ILAY,INIT,NSUB)
         INIT=0
         IF(NSUB.GT.0) WRITE(IOUT,47) PARNAM(IP)
   47    FORMAT(1X,A)
      END IF
  100 CONTINUE
C
C3------PRINT THE ARRAY.
  200 CALL ULAPRWC(ZZ,NCOL,NROW,ILAY,IOUT,IPF,ANAME)
C
C4------Return.
      RETURN
      END
      SUBROUTINE UPARARRSUB2(ZZ,NCOL,NROW,ILAY,NP,IN,IOUT,PTYP,ANAME,
     1      PACK,IPF)
C     ******************************************************************
C     Read a series of parameter names and substitute their values into
C     a 2-D array.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      double precision :: ZZ(NCOL,NROW)
      CHARACTER*(*) PTYP,PACK
      CHARACTER*24 ANAME
      CHARACTER*200 LINE
      CHARACTER*10 CTMP1,CTMP2,CTMP3,CTMP4
      double precision :: rdum
C     ------------------------------------------------------------------
C
C1------Set initialization flag to cause USUB2D to initialze ZZ to 0.
      INIT=1
C
C2------Read each parameter name.
      DO 100 N=1,NP
        READ(IN,'(A)') LINE
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
        WRITE(IOUT,5) LINE(ISTART:ISTOP)
    5   FORMAT(' Parameter:  ',A)
        IF(LINE(ISTART:ISTOP).EQ.' ') THEN
          WRITE(IOUT,*) ' Blank parameter name in the ',PACK,' file.'
          CALL USTOP(' ')
        END IF
C
C3------Loop through each parameter looking for the specified name.
        CTMP1=LINE(ISTART:ISTOP)
        CALL UPCASE(CTMP1)
        DO 10 IP=1,IPSUM
          CTMP2=PARNAM(IP)
          CALL UPCASE(CTMP2)
          IF(CTMP1.EQ.CTMP2) GO TO 20
C
C3A-----Stop looping if the end of the parameter list is found.
          IF(PARNAM(IP).EQ.' ') GO TO 15
   10   CONTINUE
   15   WRITE(IOUT,16) PACK
   16   FORMAT(1X,'Error in ',A,' file:',/
     1      1X,'The above parameter must be defined prior to its use')
        CALL USTOP(' ')
C
C4------Found parameter.
   20   CONTINUE
        IF(PARTYP(IP).NE.PTYP) THEN
C5------Print an error message if the parameter type does not match.
          WRITE(IOUT,83) PARNAM(IP),PARTYP(IP),PACK,PTYP
   83     FORMAT(1X,'Parameter type conflict:',/
     1           1X,'Named parameter:',A,' was defined as type:',A,/
     2           1X,'However, this parameter is used in the ',A,
     3             ' file, so it should be type:',A)
          CALL USTOP(' ')
        ENDIF
C
C6------Check to see if this parameter is time varying (has instances).
        NUMINST=IPLOC(3,IP)
        ILOC=IPLOC(4,IP)
        NI=1
C
C6A-----If parameter is time-varying, read instance name.
        IF(NUMINST.GT.0) THEN
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
          CTMP3=LINE(ISTART:ISTOP)
          IF(CTMP3.EQ.' ') THEN
            WRITE(IOUT,1000) PACK,PARNAM(IP)
 1000       FORMAT(/,1X,'Blank instance name in the ',A,
     &             ' file for parameter ',A)
            CALL USTOP(' ')
          ENDIF
          WRITE(IOUT,1010) CTMP3
 1010     FORMAT(3X,'Instance:  ',A)
          CALL UPCASE(CTMP3)
C
C6B------Look for instance name
          DO 50 KI=1,NUMINST
            CTMP4=INAME(ILOC+KI-1)
            CALL UPCASE(CTMP4)
            IF(CTMP3.EQ.CTMP4) THEN
              NI=KI
              GOTO 55
            ENDIF
   50     CONTINUE
          WRITE(IOUT,1020) PACK,CTMP3,PARNAM(IP)
 1020     FORMAT(/,1X,'The ',A,' file specifies undefined instance "',
     &           A,'" for parameter ',A)
          CALL USTOP(' ')
   55     CONTINUE
        ENDIF
C
C7------Check to see if this parameter is already active.
        IF (IACTIVE(IP).GT.0) THEN
          WRITE(IOUT,1030) PARNAM(IP)
 1030     FORMAT(/,1X,'*** ERROR: PARAMETER "',A,
     &        '" HAS ALREADY BEEN ACTIVATED THIS STRESS PERIOD',/,
     &        ' -- STOP EXECUTION (UPARARRSUB2)')
          CALL USTOP(' ')
        ENDIF
C
C8------Activate the parameter and substitute.  Reset INIT so that
C8------any further calls to USUB2D will not reinitilize ZZ.
        IACTIVE(IP)=NI
        II=IP
        CALL USUB2D(ZZ,NCOL,NROW,II,ILAY,INIT,NSUB)
        INIT=0
C
C9------Get new value of print flag if it is there.
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,RDUM,-1,IN)
        IF(LINE(ISTART:ISTOP) .NE.'E' .AND.
     1     LINE(ISTART:ISTOP) .NE.' ') IPF=I
C
  100 CONTINUE
C
C10-----PRINT THE ARRAY.
  200 CALL ULAPRWC(ZZ,NCOL,NROW,ILAY,IOUT,IPF,ANAME)
C
C11-----Return.
      RETURN
      END
      SUBROUTINE USUB2D(ZZ,NCOL,NROW,IP,ILAY,INIT,NSUB)
C     ******************************************************************
C     Substitute values for a single parameter into a 2-D array.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      double precision :: ZZ(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------Define constants.
      ZERO=0.0
C
C2------Initialize the array if INIT is not 0.
      IF(INIT.NE.0) THEN
        DO 10 I=1,NROW
          DO 5 J=1,NCOL
            ZZ(J,I)=ZERO
    5     CONTINUE
   10   CONTINUE
      END IF
C
C3------Identify clusters, which depends on the instance if the
C3------parameter is a time varying parameter.
      ICSTART=IPLOC(1,IP)
      ICSTOP=IPLOC(2,IP)
      NUMINST=IPLOC(3,IP)
      IF(NUMINST.GT.1) THEN
C       Select correct instance
        NCLU=(ICSTOP-ICSTART+1)/NUMINST
        NI=IACTIVE(IP)
        ICSTART=ICSTART+(NI-1)*NCLU
        ICSTOP=ICSTART+NCLU-1
      ENDIF
C
C4------Loop through each cluster definition for layers that match the
C4------specified layer.
      NSUB=0
      DO 80 IC=ICSTART,ICSTOP
C
C4A-----Check if the cluster layer matches the specified layer
        IF(IPCLST(1,IC).EQ.ILAY) THEN
C
C4B-----The parameter layer matches the specified layer.  Look at zone
C4B-----value to determine which cells to substitute. Also identify the
C4B-----multiplier array.
          MLT=IPCLST(2,IC)
          AA=1.
          IZ=IPCLST(3,IC)
          IF(IZ.GT.0) THEN
C
C4C-----IZ>0. Loop through all cells.  If the value in the zone array
C4C-----is equal to one of the cluster zone values, add the parameter
C4C-----value into the array.
            DO 50 I=1,NROW
              DO 40 J=1,NCOL
                DO 30 JJ=5,IPCLST(4,IC)
                  IF(IZON(J,I,IZ).EQ.IPCLST(JJ,IC)) THEN
                    IF(MLT.GT.0) AA=RMLT(J,I,MLT)
                    ZZ(J,I)=ZZ(J,I)+AA*B(IP)
                    NSUB=NSUB+1
                  END IF
   30           CONTINUE
   40         CONTINUE
   50       CONTINUE
          ELSE
C
C4D-----IZ is 0.  Loop through all cells adding the parameter value into
C4D-----the array.
            DO 70 I=1,NROW
              DO 60 J=1,NCOL
                IF(MLT.GT.0) AA=RMLT(J,I,MLT)
                ZZ(J,I)=ZZ(J,I)+AA*B(IP)
   60         CONTINUE
   70       CONTINUE
            NSUB=NSUB+NCOL*NROW
          END IF
        END IF
   80 CONTINUE
C
C5------Return.
      RETURN
      END
      SUBROUTINE UPARLSTAL(IN,IOUT,LINE,NP,MXL)
C     ******************************************************************
C     Setup list parameter definition for a package
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      CHARACTER*(*) LINE
      double precision :: r
C     ------------------------------------------------------------------
C
C1------Decode PARAMETER definitions if they exist
      NP=0
      MXL=0
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'PARAMETER') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NP,R,IOUT,IN)
         IF(NP.LT.0) NP=0
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MXL,R,IOUT,IN)
         IF(MXL.LT.0) MXL=0
         WRITE(IOUT,31) NP,MXL
   31    FORMAT(1X,I10,' Named Parameters     ',I10,' List entries')
         READ(IN,'(A)') LINE
      ELSE
         WRITE(IOUT,'(A)') ' No named parameters'
      END IF
C
C2------Return.
      RETURN
      END



      SUBROUTINE UPARLSTRP(LSTSUM,MXLST,IN,IOUT,NP,PACK,PTYPX,ITERP,
     &                     NUMINST,storepar)
C     ******************************************************************
C     Read and store list parameter definition information for one
C     parameter.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      logical :: storepar
      CHARACTER*(*) PACK,PTYPX
      CHARACTER*4 PTYP
      CHARACTER*10 PN,CTMP1,CTMP2
      CHARACTER*200 LINE
      double precision :: r, pv
C     ------------------------------------------------------------------
C
C1------Read the parameter name and definition.
      READ(IN,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
      PN=LINE(ISTART:ISTOP)
      CTMP1=PN
      CALL UPCASE(CTMP1)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      PTYP=LINE(ISTART:ISTOP)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,PV,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NLST,R,IOUT,IN)
C
C2------Check for multiple instances.
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF (LINE(ISTART:ISTOP).EQ.'INSTANCES') THEN
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NUMINST,R,IOUT,IN)
        IF (NUMINST.LT.1) THEN
          WRITE(IOUT,12) PN,PTYP
   12     FORMAT(/,1X,'*** ERROR: NUMINST SPECIFIED LESS THAN 1',
     &        ' FOR PARAMETER "',A,'"',/,12X,'OF TYPE "',A,
     &        '" -- STOP EXECUTION (UPARLSTRP)')
          CALL USTOP(' ')
        ENDIF
      ELSE
        NUMINST = 0
      ENDIF
C

C3------Look for parameter in list of parameters.
      if (storepar) then
      DO 10 NP=1,MXPAR
        CTMP2=PARNAM(NP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
C
C3A-----If found, determine if it is an illegal duplicate.
          IF(PARTYP(NP).NE.' ' .AND. ITERP.EQ.1) THEN
C
C3B-----Illegal duplicate.
            WRITE(IOUT,110) CTMP1
  110       FORMAT (' Duplicate parameter name: ',A)
            CALL USTOP(' ')
          END IF
C
C3C-----Parameter was predefined in SEN file (PARTYP blank) or in
C3C-----a prior simulation (ITERP not 1) -- leave its value alone.
C3C-----(i.e. ignore PV).
          GO TO 100
        ELSE IF(PARNAM(NP).EQ.' ') THEN
C
C4------Parameter was not found in the list, so it is a new definition.
          PARNAM(NP)=PN
          B(NP)=PV
          IPSUM=IPSUM+1
          GO TO 100
        ENDIF
10    CONTINUE
C
C5------Too many parameters.
      WRITE(IOUT,99)
  99  FORMAT(' Number of parameters exceeds MXPAR -- STOP EXECUTION')
      CALL USTOP(' ')
      endif
C
C6------Parameter is a new parameter, or it was predefined in the
C6------Parameter Value file or defined in a previous simulation.  Continue
C6------processing.
 100  CONTINUE
      IF(ITERP.EQ.1) THEN
C
C7------Parameter is new or was predefined in the Parameter Value file.
C7------Process the remaining parameter information.
        if (storepar) then
          PARTYP(NP)=PTYP
          IPLOC(1,NP)=LSTSUM
        endif
        NI=MAX(1,NUMINST)
        LSTSUM=LSTSUM+(NLST*NI)
        
        if (.not. storepar) return
        
        IPLOC(2,NP)=LSTSUM-1
        IPLOC(3,NP)=NUMINST
        IPLOC(4,NP)=INAMLOC
        INAMLOC=INAMLOC+NUMINST
C
C8------WRITE PARAMETER INFORMATION
        WRITE(IOUT,121) PARNAM(NP),PARTYP(NP)
  121   FORMAT(1X/,1X,'PARAMETER NAME:',A,'   TYPE:',A)
        WRITE(IOUT,122) PV
  122   FORMAT(1X,'Parameter value from package file is: ',1PG13.5)
        IF(B(NP).NE.PV) THEN
          WRITE(IOUT,123) B(NP)
  123     FORMAT(1X,'This value has been changed to:',7X,1PG13.5,
     &        ', as read from',/,' the Parameter Value file')
        END IF
        WRITE(IOUT,130) NLST
  130   FORMAT(  '   NUMBER OF ENTRIES: ',I6)
        IF(NUMINST.GT.0) THEN
          WRITE(IOUT,131)NUMINST
  131     FORMAT('   NUMBER OF INSTANCES: ',I4)
        ENDIF
C
C9------Check if the parameter list will fit in the package list array.
        IF((LSTSUM-1) .GT. MXLST) THEN
          WRITE(IOUT,134) LSTSUM-1,MXLST
  134     FORMAT(1X,'EXCEEDED THE MAXIMUM NUMBER OF LIST ENTRIES:'/
     1          1X,I5,' list entries have been specified'/
     2          1X,'The maximum number of list entries is',I5)
          CALL USTOP(' ')
        END IF
C
C10-----Check if number of instances exceeds the maximum allowed.
        IF((INAMLOC-1).GT.MXINST) THEN
          WRITE(IOUT,135)INAMLOC-1,MXINST
  135     FORMAT(1X,'EXCEEDED THE MAXIMUM NUMBER OF INSTANCES:'/
     1          1X,I5,' instances have been specified'/
     2          1X,'The maximum number of instances is',I5)
          CALL USTOP(' ')
        END IF
C
C11-----Check for correct parameter type.
        IF(PARTYP(NP).NE.PTYPX) THEN
          WRITE(IOUT,137) PTYPX,PACK
  137     FORMAT(1X,'Parameter type must be:',A,' in the ',A,' Package')
          CALL USTOP(' ')
        END IF
C
C12-----Parameter definition must include at least one cell.
        IF (NLST.LE.0) THEN
          WRITE(IOUT,140) PN
          CALL USTOP(' ')
        ENDIF
  140   FORMAT(' ERROR:  DEFINITION FOR PARAMETER "',A,'"',
     &      ' INCLUDES NO CELLS',/,'   -- STOP EXECUTION (UPARLSTRP)')
      ELSE
C
C13-----This is not the first time the simulation was run, so the parameter
C13-----was already defined.  Set values of arguments to be returned.
        LSTSUM=LSTSUM+IPLOC(2,NP)-IPLOC(1,NP)+1
        NUMINST=IPLOC(3,NP)
      ENDIF
C
C14-----Set the parameter to be inactive.
      IACTIVE(NP)=0
C
C15-----Return.
      RETURN
      END



      SUBROUTINE UINSRP(I,IN,IOUT,IP,ITERP)
C     ******************************************************************
C     Read and store one instance name.
C     I is the instance number, and IP is the parameter number.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      CHARACTER*10 CTMP1,CTMP2
      CHARACTER*200 LINE
      double precision :: r
C     ------------------------------------------------------------------
C
C1------COMPUTE LOCATION OF NAME IN INAME, AND READ LINE CONTAINING
C1------INSTANCE NAME.
      IPL4 = IPLOC(4,IP)
      ILOC = IPL4+I-1
      READ(IN,1000) LINE
 1000 FORMAT(A)
C
C2------GET INSTANCE NAME AND STORE IN INAME.
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
      INAME(ILOC) = LINE(ISTART:ISTOP)
      CTMP1 = LINE(ISTART:ISTOP)
      CALL UPCASE(CTMP1)
C
C3------WRITE NAME UNLESS THIS IS NOT THE FIRST TIME THE SIMULATION
C3------HAS BEEN RUN.
      IF(ITERP.EQ.1) WRITE(IOUT,1010)INAME(ILOC)
 1010 FORMAT(/,3X,'INSTANCE:  ',A)
C
C4------CHECK FOR DUPLICATE INSTANCE NAME IF THIS IS NOT THE FIRST
C4------INSTANCE.
      IF (I.GT.1) THEN
        DO 10 J=IPL4,IPL4+I-2
          CTMP2 = INAME(J)
          CALL UPCASE(CTMP2)
          IF (CTMP1.EQ.CTMP2) THEN
            WRITE(IOUT,1020)INAME(J)
 1020       FORMAT(/,1X,'*** ERROR: "',A,
     &'" IS A DUPLICATE INSTANCE NAME FOR THIS PARAMETER',/,
     &' -- STOP EXECUTION (UINSRP)')
            CALL USTOP(' ')
          ENDIF
   10   CONTINUE
      ENDIF
C
C5------RETURN.
      RETURN
      END
      SUBROUTINE UPARLSTSUB(IN,PACK,IOUTU,PTYP,RLIST,LSTVL,LSTDIM,NREAD,
     1                MXLST,NTOT,IPVL1,IPVL2,LABEL,CAUX,NCAUX,NAUX)
C     ******************************************************************
C     Read a list parameter name, look it up in the list of parameters,
C     and substitute values into active part of package array.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      CHARACTER*(*) PACK,PTYP
      double precision :: RLIST(LSTVL,LSTDIM)
      CHARACTER*(*) LABEL
      CHARACTER*16 CAUX(NCAUX)
      CHARACTER*200 LINE
      CHARACTER*10 CTMP1,CTMP2,CTMP3,CTMP4
      double precision :: rdum
C     ------------------------------------------------------------------
C
C1------The LIST file unit is the absolute value of IOUTU.
C1------Read the parameter name.
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
C
C2------Find the parameter in the list of parameters.
      CTMP1=LINE(ISTART:ISTOP)
      CALL UPCASE(CTMP1)
      DO 100 IP=1,IPSUM
        CTMP2=PARNAM(IP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
          IF(PARTYP(IP).NE.PTYP) THEN
            WRITE(IOUT,11) PARNAM(IP),PARTYP(IP),PACK,PTYP
   11       FORMAT(1X,'Parameter type conflict:',/
     1        1X,'Named parameter:',A,' was defined as type:',A,/
     2        1X,'However, this parameter is used in the ',A,
     3          ' file, so it should be type:',A)
            CALL USTOP(' ')
          END IF
C
C3------Set indices to point to the cells that correspond to the
C3------specified parameter.  If the parameter is time varying, set the
C3------indices to the specified instance.
          NLST=IPLOC(2,IP)-IPLOC(1,IP)+1
          NUMINST=IPLOC(3,IP)
          ILOC=IPLOC(4,IP)
          NI=1
          IF(NUMINST.GT.0) THEN
            NLST=NLST/NUMINST
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
            CTMP3=LINE(ISTART:ISTOP)
            IF(CTMP3.EQ.' ') THEN
              WRITE(IOUT,15)PACK,PARNAM(IP)
   15         FORMAT(/,1X,'Blank instance name in the ',A,
     &               ' file for parameter ',A)
              CALL USTOP(' ')
            ENDIF
            WRITE(IOUT,17) CTMP3
   17       FORMAT(3X,'Instance:  ',A)
            CALL UPCASE(CTMP3)
            DO 50 KI=1,NUMINST
              CTMP4=INAME(ILOC+KI-1)
              CALL UPCASE(CTMP4)
              IF(CTMP3.EQ.CTMP4) THEN
                NI=KI
                GOTO 55
              ENDIF
   50       CONTINUE
            WRITE(IOUT,53) PACK,CTMP3,PARNAM(IP)
   53       FORMAT(/,1X,'The ',A,' file specifies undefined instance "',
     &             A,'" for parameter ',A)
            CALL USTOP(' ')
   55       CONTINUE
          ENDIF
C
C4------Check that the parameter is not already active.
          IF (IACTIVE(IP).GT.0) THEN
            WRITE(IOUT,73) PARNAM(IP)
   73       FORMAT(/,1X,'*** ERROR: PARAMETER "',A,
     &          '" HAS ALREADY BEEN ACTIVATED THIS STRESS PERIOD',/,
     &          ' -- STOP EXECUTION (UPARLSTSUB)')
            CALL USTOP(' ')
          ENDIF
C
C5------Set the active flag.
          IACTIVE(IP)=NI
C
C6------Accumulate the total number of active cells in the list.
          NTOT=NTOT+NLST
          IF(NTOT.GT.MXLST) THEN
            WRITE(IOUT,83) NTOT,MXLST
   83       FORMAT(1X,/1X,'THE NUMBER OF ACTIVE LIST ENTRIES (',I6,
     1       ')',/1X,'IS GREATER THAN THE MAXIMUM ALLOWED (',I6,')')
            CALL USTOP(' ')
          END IF
C
C7------Write label for list values if IOUTU is positive.
          IF (IOUTU.GT.0) CALL ULSTLB(IOUT,LABEL,CAUX,NCAUX,NAUX)
C
C8------Copy the values from the paramter location into the front part
C8------of the list where the currently active list is kept.
          DO 90 I=1,NLST
            II=NTOT-NLST+I
            III=I-1+IPLOC(1,IP)+(NI-1)*NLST
            DO 85 J=1,NREAD
              RLIST(J,II)=RLIST(J,III)
   85       CONTINUE
C
C8A-----Scale the RLIST values from IPVL1 to IPVL2 by the parameter
C8A-----value.
            DO 86 IPVL=IPVL1,IPVL2
              RLIST(IPVL,II)=RLIST(IPVL,II)*B(IP)
   86       CONTINUE
            IL=RLIST(1,II)
            IR=RLIST(2,II)
            IC=RLIST(3,II)
            IF (IOUTU.GT.0) WRITE(IOUT,89) II,IL,IR,IC,
     &          (RLIST(JJ,II),JJ=4,NREAD)
   89       FORMAT(1X,I6,I7,I7,I7,14G16.4)
   90     CONTINUE
C
C8B------After moving the data, return.
          RETURN
        END IF
  100 CONTINUE
C
C9------All parameter names have been checked without finding the
C9------parameter. Write an error message and stop.
      WRITE(IOUT,*) ' The ',PACK,
     1   ' file specifies an undefined parameter:',LINE(ISTART:ISTOP)
      CALL USTOP(' ')
C
      END
      SUBROUTINE PRESET(PTYP)
C     ******************************************************************
C     Clear active flag for all parameters of a specified type
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      CHARACTER*(*) PTYP
C     ------------------------------------------------------------------
C
C1------Loop through all parameters.  Set IACTIVE to 0 when the
C1------parameter type matches.
      DO 10 I=1,IPSUM
      IF(PARTYP(I).EQ.PTYP) IACTIVE(I)=0
   10 CONTINUE
C
C2------Return.
      RETURN
      END
      SUBROUTINE UPARLSTLOC(IN,PACK,IOUT,PTYP,IBEG,IEND,PV)
C     ******************************************************************
C     Read a list parameter name, look it up in the list of
C     parameters, and return its location.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      CHARACTER*(*) PACK,PTYP
      CHARACTER*200 LINE
      CHARACTER*10 CTMP1,CTMP2,CTMP3,CTMP4
      double precision :: rdum
C     ------------------------------------------------------------------
C
C1------Read the parameter name.
      READ(IN,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
      WRITE(IOUT,500) LINE(ISTART:ISTOP)
  500 FORMAT(' Parameter:  ',A)
      IF(LINE(ISTART:ISTOP).EQ.' ') THEN
        WRITE(IOUT,*) ' Blank parameter name in the ',PACK,' file.'
        CALL USTOP(' ')
      END IF
C
C2------Find the parameter in the list of parameters.
      CTMP1=LINE(ISTART:ISTOP)
      CALL UPCASE(CTMP1)
      DO 100 IP=1,IPSUM
        CTMP2=PARNAM(IP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
          IF(PARTYP(IP).NE.PTYP) THEN
            WRITE(IOUT,510) PARNAM(IP),PARTYP(IP),PACK,PTYP
  510       FORMAT(1X,'Parameter type conflict:',/
     &       1X,'Named parameter:',A,' was defined as type:',A,/
     &       1X,'However, this parameter is used in the ',A,
     &          ' file, so it should be type:',A)
            CALL USTOP(' ')
          END IF
C
C3------Set indices to point to the cells that correspond to the
C3------specified parameter.  If the parameter is time varying, set the
C3------indices to the specified instance.
          NLST=IPLOC(2,IP)-IPLOC(1,IP)+1
          NUMINST=IPLOC(3,IP)
          ILOC=IPLOC(4,IP)
          NI=1
          IF(NUMINST.GT.0) THEN
            NLST=NLST/NUMINST
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,IDUM,RDUM,IOUT,IN)
            CTMP3=LINE(ISTART:ISTOP)
            IF(CTMP3.EQ.' ') THEN
              WRITE(IOUT,520)PACK,PARNAM(IP)
  520         FORMAT(/,1X,'Blank instance name in the ',A,
     &       ' file for parameter ',A)
              CALL USTOP(' ')
            ENDIF
            WRITE(IOUT,530) CTMP3
  530       FORMAT(3X,'Instance:  ',A)
            CALL UPCASE(CTMP3)
            DO 50 KI=1,NUMINST
              CTMP4=INAME(ILOC+KI-1)
              CALL UPCASE(CTMP4)
              IF(CTMP3.EQ.CTMP4) THEN
                NI=KI
                GOTO 60
              ENDIF
   50       CONTINUE
            WRITE(IOUT,540) PACK,CTMP3,PARNAM(IP)
  540       FORMAT(/,1X,'The ',A,' file specifies undefined instance "',
     &       A,'" for parameter ',A)
            CALL USTOP(' ')
   60       CONTINUE
          ENDIF
C
C4------Check that the parameter is not already active..
          IF (IACTIVE(IP).GT.0) THEN
            WRITE(IOUT,550) PARNAM(IP)
  550       FORMAT(/,1X,'*** ERROR: PARAMETER "',A,
     &'" HAS ALREADY BEEN ACTIVATED THIS STRESS PERIOD',/,
     &' -- STOP EXECUTION (UPARLSTLOC)')
            CALL USTOP(' ')
          ENDIF
C
C5------Set the active flag.
          IACTIVE(IP)=NI
C
C6------Compute pointers to the beginning and ending location of the
C6------parameter list within the package list of data.
          IBEG=IPLOC(1,IP)+(NI-1)*NLST
          IEND=IBEG+NLST-1
          PV=B(IP)
          IACTIVE(IP)=NI
C
C7------Return.
          RETURN
        END IF
  100 CONTINUE
C
C8------All parameter names have been checked without finding the
C8------parameter. Write an error message and stop.
      WRITE(IOUT,*) ' The ',PACK,
     &   ' file specifies an undefined parameter:',LINE(ISTART:ISTOP)
      CALL USTOP(' ')
C
      END
      SUBROUTINE UPARARRCK(BUFF,IBOUND,IOUT,LAY,NCOL,NLAY,NROW,PTYP)
C     ******************************************************************
C     CHECK FOR COMPLETE DEFINITION OF ONE LAYER OF CELLS BY ARRAY
C     PARAMETERS OF A GIVEN TYPE.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      INTEGER IBOUND(NCOL,NROW,NLAY)
      double precision :: BUFF(NCOL,NROW)
      CHARACTER*4 PTYP
C     ------------------------------------------------------------------
C
C1------Make sure that the parameter type is non-blank.
      IF (PTYP.EQ.' ') THEN
        WRITE (IOUT,500)
  500   FORMAT(1X,'ERROR: BLANK PARAMETER TYPE -- STOP EXECUTION',
     &       ' (UPARARRCK)')
        CALL USTOP(' ')
      ENDIF
C
C2------Initialize BUFF to 0.
      DO 20 I = 1, NROW
        DO 10 J = 1, NCOL
          BUFF(J,I) = 0.0
   10   CONTINUE
   20 CONTINUE
C
C3------Loop through parameters to find matching parameter type.
C3------Increment BUFF for each cell where a parameter of the specified
C3------type applies.
      DO 100 IP = 1, IPSUM
        IF (PARTYP(IP).EQ.PTYP) THEN
C
C3A-----Loop through clusters associated with this parameter.
          DO 80 IC = IPLOC(1,IP), IPLOC(2,IP)
            IF (IPCLST(1,IC).EQ.LAY) THEN
              IZA = IPCLST(3,IC)
              DO 60 I = 1, NROW
                DO 50 J = 1,NCOL
                  IF (IZA.GT.0) THEN
C
C3B-----Loop through zones listed for this cluster.
                    DO 40 IZI = 5, IPCLST(4,IC)
                      IZ = IPCLST(IZI,IC)
                      IF (IZ.EQ.IZON(J,I,IZA)) THEN
                        BUFF(J,I) = BUFF(J,I) + 1.0
                      ENDIF
   40               CONTINUE
                  ELSE
C
C3C-----Zones do not apply to this cluster -- apply to all cells.
                    BUFF(J,I) = BUFF(J,I) + 1.0
                  ENDIF
   50           CONTINUE
   60         CONTINUE
            ENDIF
   80     CONTINUE
        ENDIF
  100 CONTINUE
C
C4------Identify any active cells where BUFF is equal to zero, which
C4------indicates cells that are not defined by any parameter of the
C4------specified type applies.
      IERR = 0
      DO 140 I = 1, NROW
        DO 120 J = 1, NCOL
          IF (IBOUND(J,I,LAY).NE.0) THEN
            IF (BUFF(J,I).EQ.0.0)THEN
              WRITE (IOUT,510) I,J,LAY,PTYP
  510         FORMAT(1X,'ROW: ',I5,', COLUMN: ',I5,' IN LAYER ',I3,
     &       ' NOT DEFINED FOR PARAMETER TYPE ',A)
              IERR = IERR + 1
            ENDIF
          ENDIF
  120   CONTINUE
  140 CONTINUE
C
C5------IF any active cells were found with undefined values, write an
C5------error message and stop.
      IF (IERR.GT.0) THEN
        WRITE (IOUT,520)
  520   FORMAT(/,1X,'PARAMETER DEFINITIONS INCOMPLETE -- STOP',
     &       ' EXECUTION (UPARARRCK)')
        CALL USTOP(' ')
      ENDIF
C
C6------Return.
      RETURN
      END
      SUBROUTINE UPARFIND(PNAME,PTYP,CPACK,IFOUND,IOUT)
C     ******************************************************************
C     Find the parameter number for a parameter name
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE PARAMMODULE
      use utl7module
      use SimPHMFModule
      CHARACTER*(*) PNAME,PTYP,CPACK
      CHARACTER*10 CTMP1,CTMP2
C     ------------------------------------------------------------------
C
C1------Abort if parameter name is blank.
      IF(PNAME.EQ.' ') THEN
         WRITE(IOUT,*) ' Blank parameter name in the ',CPACK,' file.'
         CALL USTOP(' ')
      END IF
C
C2------Look for the parameter name in the list of parameters.
      CTMP1=PNAME
      CALL UPCASE(CTMP1)
      DO IP=1,MXPAR
        CTMP2=PARNAM(IP)
        CALL UPCASE(CTMP2)
        IF(CTMP1.EQ.CTMP2) THEN
           IF(PARTYP(IP).NE.PTYP) THEN
              WRITE(IOUT,10) PARNAM(IP),PARTYP(IP),CPACK,PTYP
              CALL USTOP(' ')
           END IF
           IFOUND=IP
           RETURN
        END IF
      END DO
C
C3------Failed to find the name -- abort.
      WRITE(IOUT,20) CPACK
   10 FORMAT(1X,'Parameter type conflict:',/,1X,'Named parameter:',A,
     1       ' was defined as type:',A,/,1X,
     2       'However, this parameter is used in the ',A,
     3       ' file, so it should be type:',A)
   20 FORMAT(1X,'Parameter for ',A,' Package has not been defined')
      CALL USTOP(' ')
C
      END
