      module GwfLgrSubsModule
        
        use ConstantsModule, only: MAXCHARLEN
        use GlobalVariablesModule, only: LgrBilinear, optfile
        use InputOutputModule, only: GetUnit
        use LGRMODULE, only: SGWF2LGR2PNT, SGWF2LGR2PSV
        use OpenSpecModule, only: ACCESS, ACTION, FORM
        use SimPHMFModule, only: ustop
        use UtilitiesModule, only: GetArgs
        use utl7module, only: URDCOM, URWORD

      contains
      
      SUBROUTINE GETNAMFIL(ILGR,NGRIDS,FNAME,ILUNIT,basnam)
C     ******************************************************************
C     GET THE NAME OF THE NAME FILE
C     ******************************************************************
C        SPECIFICATIONS:
C
C     ------------------------------------------------------------------
      implicit none
      integer :: ilgr, ngrids, ilunit, icol, istart, istop, n, nc
      CHARACTER*(*) FNAME
      CHARACTER*200 COMLIN,LINE
      character(len=MAXCHARLEN) :: namfil, basnam
      LOGICAL EXISTS
      double precision :: r
C     ------------------------------------------------------------------
C
C Get name file from command line or user interaction.
        FNAME=' '
        COMLIN=' '
        basnam = ''
        call GetArgs(namfil, basnam)
        ICOL = 1
        IF(namfil.NE.' ') THEN
          FNAME = namfil
        ELSE
   15     continue
          WRITE (*,*) ' Enter the name of the NAME FILE or LGR ',
     &                'CONTROL FILE:'
          READ (*,'(A)') FNAME
          CALL URWORD(FNAME,ICOL,ISTART,ISTOP,0,N,R,0,0)
          FNAME=FNAME(ISTART:ISTOP)
          IF (FNAME.EQ.' ') GOTO 15
!          write(*,*)' Enter the name of an MF5to6 Options file:'
!          read(*,'(a)') optfile
        ENDIF
        INQUIRE (FILE=FNAME,EXIST=EXISTS)
        IF(.NOT.EXISTS) THEN
          NC=INDEX(FNAME,' ')
          FNAME(NC:NC+3)='.nam'
          INQUIRE (FILE=FNAME,EXIST=EXISTS)
          IF(.NOT.EXISTS) THEN
            WRITE (*,480) FNAME(1:NC-1),FNAME(1:NC+3)
  480       FORMAT(1X,'Can''t find name file ',A,' or ',A)
            CALL USTOP(' ')
          ENDIF
        ENDIF
C1A-----CHECK FOR LGR KEYWORD.  IF AN LGR SIMULATION, THEN READ THE 
C1A-----NUMBER OF GRIDS AND LEAVE FILE OPEN FOR PARSING.
C1A-----IF NOT LGR, THEN CLOSE FILE AND CONTINUE AS NORMAL.
      ilunit = GetUnit()
      OPEN(UNIT=ILUNIT,FILE=FNAME,STATUS='OLD',ACTION=ACTION(1))
      CALL URDCOM(ILUNIT,0,LINE)
      ICOL=1
      CALL URWORD(LINE,ICOL,ISTART,ISTOP,1,N,R,0,ILUNIT)
      IF(LINE(ISTART:ISTOP) .EQ. 'LGR') THEN
        ILGR = 1
        WRITE(*,*) ' LGR is active '
        call URWORD(line,icol,istart,istop,1,n,r,0,ilunit)
        if (line(istart:istop) == 'TRIANGULAR') then
          LgrBilinear = .false.
        elseif (line(istart:istop) == 'BILINEAR') then
          LgrBilinear = .true.
        endif
        CALL URDCOM(ILUNIT,0,LINE)
        ICOL=1
        CALL URWORD(LINE,ICOL,ISTART,ISTOP,2,NGRIDS,R,0,ILUNIT)
        WRITE(*,*) 'NGRIDS = ', NGRIDS 
      ELSE
        ilgr = 0
        CLOSE(ILUNIT)
        ilunit = 0
      ENDIF
C
      RETURN
      END SUBROUTINE GETNAMFIL

C***********************************************************************

C-----VERSION 1.0 15FEBRUARY2006 GETNAMFILLGR
      SUBROUTINE GETNAMFILLGR(INLGR,FNAME,IGRID)
C     ******************************************************************
C     READ NAMES OF THE CORRESPONDING NAME FILES FROM LGR CONTROL FILE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*200 LINE, FNAME
      LOGICAL EXISTS
      double precision :: r
C     ------------------------------------------------------------------
C1-----READ IN THE NAME OF THE NAME FILE FOR THIS GRID
      CALL URDCOM(INLGR,0,LINE)
      ICOL = 1
      CALL URWORD(LINE,ICOL,ISTART,ISTOP,0,N,R,0,0)
      FNAME=LINE(ISTART:ISTOP)
      INQUIRE (FILE=FNAME,EXIST=EXISTS)
      IF(.NOT.EXISTS) THEN
        NC=INDEX(FNAME,' ')
        FNAME(NC:NC+3)='.nam'
        INQUIRE (FILE=FNAME,EXIST=EXISTS)
        IF(.NOT.EXISTS) THEN
          WRITE (*,480) FNAME(1:NC-1),FNAME(1:NC+3)
  480     FORMAT(1X,'Can''t find name file ',A,' or ',A)
          CALL USTOP(' ')
        ENDIF
      ENDIF
C
      RETURN
      END SUBROUTINE GETNAMFILLGR
      
C***********************************************************************

C-----VERSION 1.1 10OCTOBER2006 GWF2LGR1AR      
      SUBROUTINE GWF2LGR2AR(INLGR,FNAME,NGRIDS,IGRID)
C     ******************************************************************
C     ALLOCATE SPACE FOR LOCAL GRID REFINEMENT AND READ DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IOUT,IUNIT,IBOUND,GLOBALDAT
      USE LGRMODULE,   ONLY:ISCHILD,NGRDS,NPLBEG,NPRBEG,NPCBEG,NPLEND,
     1                      NPREND,NPCEND,NCPP,NPL,IBOTFLG,ISHFLG,IBFLG,
     2                      IUPBHSV,IUCBHSV,IUPBFSV,IUCBFSV,MXLGRITER,
     3                      IOUTLGR,NBNODES,NPBNODES,IBMAXH,IBMAXF,
     4                      NCMAXH,NCMAXF,RELAXH,RELAXF,HCLOSELGR,
     5                      FCLOSELGR,HDIFFM,FDIFFM,PRATIN,CRATIN,
     6                      PRATOUT,CRATOUT,IBPFLG,IEDG,JEDG,NCPPL,
     7                      NODEH,NODEF,NCON,KPLC,IPLC,JPLC,IFACEGN,
     8                      ICBOUND,GNHEAD,DHGN,GNFLUX,GNFLUXR,
     9                      GNFLUXOLD,HOLDC,GNCOND,VCB,HK,
     &                      VK,LGRDAT 
C
      CHARACTER*200 LINE, FNAME 
      CHARACTER*14 GRIDSTATUS
      double precision :: r
C     ------------------------------------------------------------------
      ALLOCATE(ISCHILD,NGRDS,NPLBEG,NPRBEG,NPCBEG,NPLEND,NPREND,NPCEND,
     1         NCPP,NPL,IBOTFLG,ISHFLG,IBFLG,IUPBHSV,IUCBHSV,IUPBFSV,
     2         IUCBFSV,MXLGRITER,IOUTLGR,NBNODES,NPBNODES,IBMAXH,IBMAXF,
     3         NCMAXH,NCMAXF,RELAXH,RELAXF,HCLOSELGR,FCLOSELGR,HDIFFM,
     4         FDIFFM,IBPFLG(NGRIDS))
C1------PRINT A MESSAGE IDENTIFYING LGR PACKAGE
      WRITE(IOUT,500)IGRID,TRIM(FNAME)
  500 FORMAT(1X,/1X,'LGR2 -- LOCAL GRID REFINEMENT',
     &    ', VERSION 2.0, 06/25/2013',/8X,'INPUT READ FOR MODEL ',
     &    I2,' DEFINED BY NAME FILE ',A)

      ZERO = 0.
C
C2------INITIALIZE IBPLG TO ZERO
      DO LG=1,NGRIDS
        IBPFLG(LG) = 0
      ENDDO
C
C3------READ IN GRIDSTATUS, SET ISCHILD, AND PRINT COMMENTS
      CALL URDCOM(INLGR,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INLGR)
      GRIDSTATUS = LINE(ISTART:ISTOP)
      NGRDS=NGRIDS
C
C3A-----PRINT COMMENTS
      IF(GRIDSTATUS .EQ. 'PARENTONLY')THEN
        ISCHILD = -1
        WRITE(IOUT,410) 
      ELSEIF(GRIDSTATUS .EQ. 'CHILDONLY')THEN
        ISCHILD = 1
        WRITE(IOUT,411) 
      ELSEIF(GRIDSTATUS .EQ. 'PARENTANDCHILD')THEN
        ISCHILD = 0
        WRITE(IOUT,412) 
      ELSE
        WRITE(IOUT,413) GRIDSTATUS
        STOP
      ENDIF
  410 FORMAT (1X,/1X,'LOCAL GRID REFINEMENT IS ACTIVE FOR PARENT ONLY')
  411 FORMAT (1X,/1X,'LOCAL GRID REFINEMENT IS ACTIVE FOR CHILD ONLY')
  412 FORMAT (1X,'LOCAL GRID REFINEMENT IS ACTIVE FOR PARENT AND',
     &        ' CHILD') 
  413 FORMAT (1X,'INVALID INPUT FOR GRIDSTATUS IN LGR INPUT FILE: ',A,/,
     &        1X,'GRIDSTATUS MUST BE PARENTONLY, CHILDONLY, OR ',
     &         'PARENTANDCHILD')
C
C4------READ, PRINT, AND CHECK LGR DATA
C4A-----IF A PARENT GRID, READ IN OPTIONS FOR SAVING BOUNDARY HEADS 
C4A-----AND FLUXES
      IF(ISCHILD .LT. 0)THEN
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  IUPBHSV,IUPBFSV
        IF(IUPBHSV .NE. 0) WRITE(IOUT,512) IUPBHSV
        IF(IUPBFSV .NE. 0) WRITE(IOUT,513) IUPBFSV
      ENDIF
C4B-----CHECK IF A CHILD GRID
C4B-----READ IN ISHFLG,IUCBHSV,IUCBFSV,MXLGRITER,RELAXH,RELAXF,HCLOSELGR,
C4B-----AND FCLOSELGR.  STORE IBFLG IN IBPFLG OF THE PARENT GRID.
      IF(ISCHILD .GE. 0)THEN
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  ISHFLG,IBFLG,IUCBHSV,IUCBFSV 
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  MXLGRITER,IOUTLGR
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  RELAXH,RELAXF
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  HCLOSELGR,FCLOSELGR
C  NOTE:  THIS IS HARDWIRED FOR A SINGLE PARENT GRID.  THE INDEX SHOULD
C  BE INDEXED TO WHATEVER THE CURRENT PARENT IS FOR THIS CHILD.
        LGRDAT(1)%IBPFLG(IGRID)=IBFLG

C4B1----PRINT COMMENTS
        IF(ISHFLG .EQ. 1) WRITE(IOUT,510) ISHFLG
        IF(IBFLG  .GT. 0)THEN
          WRITE(IOUT,511) IBFLG
        ELSE
          CALL USTOP('IBFLG MUST BE > 0 FOR CHILD GRID')
        ENDIF
        IF(IUCBHSV .NE. 0) WRITE(IOUT,512) IUCBHSV
        IF(IUCBFSV .NE. 0) WRITE(IOUT,513) IUCBFSV
        WRITE(IOUT,515) MXLGRITER
        IF(IOUTLGR .LT. 0) WRITE(IOUT,516) IOUTLGR
        IF(IOUTLGR .EQ. 0) WRITE(IOUT,517) IOUTLGR
        IF(IOUTLGR .GT. 0) WRITE(IOUT,518) IOUTLGR
        WRITE(IOUT,520) RELAXH,RELAXF
        WRITE(IOUT,521) HCLOSELGR,FCLOSELGR
cswm        IF(RELAXH .LE. ZERO .OR. RELAXF .LE. ZERO)THEN
        IF(RELAXH .LE. ZERO)THEN
          WRITE(IOUT,*) 'RELAXATION FACTOR RELAXH MUST BE > 0'
          CALL USTOP(' ')
        ENDIF
        IF(RELAXF .LE. ZERO) RELAXF =1. ! kludge to allow new formulation
  510   FORMAT (15X, 'STARTING HEADS FROM PARENT WILL BE USED:',
     &              ' ISHFLG = ',I3)
  511   FORMAT (15X,'VALUE IN IBOUND INDICATING BOUNDARY ',
     &                  'INTERFACE = ', I4)
  512   FORMAT (15X,'BOUNDARY HEADS WILL BE SAVED ON UNIT ', I3)
  513   FORMAT (15X,'BOUNDARY FLUXES WILL BE SAVED ON UNIT ', I3)
  515   FORMAT (15X,'MAX NUMBER OF LGR ITERATIONS =', I3)
  516   FORMAT (15X,'LGR ITERATION RESULTS WRITTEN TO SCREEN: ',
     &                   'IOUTLGR = ', I3)
  517   FORMAT (15X,'LGR ITERATION RESULTS NOT WRITTEN: ',
     &                   'IOUTLGR = ', I3)
  518   FORMAT (15X,'LGR ITERATION RESULTS WRITTEN TO FILE: ',
     &                   'IOUTLGR = ', I3)
  520   FORMAT (1X,/15X,'WEIGHTING FACTORS FOR RELAXATION',/,
     &            15x,'RELAXH(HEAD)' ,1x,'RELAXF(FLUX)',/,15X,25('-'),/,
     &            14x,E10.3,5X,E10.3) 
  521   FORMAT (1X,/15x,'CLOSURE CRITERIA FOR LGR ITERATIONS',/,
     &            15x,'HCLOSELGR' ,7x,'FCLOSELGR',/,15X,25('-'),/,14x,
     &            1P,E10.3,5X,E10.3)
C
C4C-----READ IN DATA FOR CHILD/PARENT MESH INTERFACE LOCATION AND 
C4C-----REFINEMENT INFO
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  NPLBEG,NPRBEG,NPCBEG
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  NPLEND,NPREND,NPCEND  
        CALL URDCOM(INLGR,IOUT,LINE)
        READ (LINE,*)  NCPP
C4C1----ALLOCATE AND READ LAYER REFINEMENT
        NPL = NPLEND - NPLBEG + 1
        ALLOCATE(NCPPL(NPL))
        READ(INLGR,*) (NCPPL(K),K=1,NPL)
C
C4D-----PRINT REFINEMENT INFO
        WRITE(IOUT, 525) NPLBEG,NPRBEG,NPCBEG,
     &                   NPLEND,NPREND,NPCEND
  525   FORMAT (1x,/,15x,'STARTING LAYER, ROW, COLUMN=',I4,',',I4,',',
     &        I4,/15x,'ENDING LAYER, ROW, COLUMN=  ',I4,',',I4,',',I4)
        WRITE(IOUT,530) NCPP
  530   FORMAT (1x,14x,'NCPP: NUMBER OF CHILD CELLS PER WIDTH OF',
     &                 ' PARENT CELL=',I2)
        WRITE(IOUT,535) (K,NCPPL(K),K=1,NPL)
  535   FORMAT (1x,14X,'NCPPL: NUMBER OF CHILD LAYERS IN LAYER ',I2,
     &                 ' OF PARENT =',I2)
C
C4D1----CHECK REFINEMENT RATIOS AND GRID DISCRETIZATION
        IF(NPLBEG .NE. 1)THEN
          WRITE(IOUT,*)'NPLBEG IS NOT = 1  REFINEMENT MUST BEGIN IN ',
     &                 'TOP LAYER'
          CALL USTOP(' ')
        ENDIF
        IF(NCPP*(NPREND-NPRBEG + 1) .NE. NROW)THEN
          WRITE(IOUT,*)'NROW DOES NOT EQUAL NCPP*(NPREND - NPRBEG + 1)'
          CALL USTOP(' ')
        ELSEIF(NCPP*(NPCEND-NPCBEG + 1) .NE. NCOL)THEN
          WRITE(IOUT,*)'NCOL DOES NOT EQUAL NCPP*(NPCEND - NPCBEG + 1)'
          CALL USTOP(' ')
        ENDIF
C
C4D2----CHECK VERTICAL REFINEMENT 
C4D2----HARDWIRED TO TOP PARENT. IT SHOULD BE INDEXED TO 
C4D2----CURRENT PARENT FOR THIS GRID.
        IBOTFLG = 1
        IF(NPLEND .EQ. GLOBALDAT(1)%NLAY) IBOTFLG=0
        NSUM = 0
        DO K=1,NPL
          NSUM = NSUM + NCPPL(K) 
        ENDDO
        IF(NSUM .NE. NLAY)THEN
          WRITE(IOUT,*)'VERTICAL REFINEMENT DOES NOT ALIGN', 
     &                 ' WITH NLAY: CHECK NCPPL'
          CALL USTOP(' ')
        ENDIF

C END CHILD GRID DATA
      ENDIF
C

C     
C  DETERMINE IF THIS IS A PARENT, CHILD, OR CHILD/PARENT GRID 
C  AND ALLOCATE AS SUCH.
      IF(ISCHILD .GE. 0)THEN
        ALLOCATE(PRATIN,CRATIN,PRATOUT,CRATOUT)
        ALLOCATE(NODEH(3),NODEF(3),ICBOUND(NCOL,NROW,NLAY))
C SET NUMBER OF BOUNDARY NODES FOR STORING CHILD TO PARENT LOCATIONS
        ICBOUND=ABS(IBOUND)
        NBNODES = COUNT(ICBOUND .EQ. IBFLG) 
C
C  ALLOCATE SPACE FOR MAPPING CHILD LOCATIONS TO PARENT LOCATIONS
        ALLOCATE(NCON(NBNODES),KPLC(3,NBNODES),IPLC(3,NBNODES),
     &           JPLC(3,NBNODES),IFACEGN(3,NBNODES),IEDG(NBNODES),
     &           JEDG(NBNODES))
        ALLOCATE(HOLDC(NCOL,NROW,NLAY),
     &           GNHEAD(3,NBNODES),DHGN(3,NBNODES),GNFLUX(3,NBNODES),
     &           GNFLUXR(3,NBNODES),GNFLUXOLD(3,NBNODES),
     &           GNCOND(3,NBNODES))
C
        ALLOCATE(VCB(4))
      ENDIF
C
C-------ALLOCATE SPACE FOR STORING HK AND VK DEPENDING ON FLOW PACKAGE
C-------NOTE: ALLOCATING NLAY FOR VK EVEN THOUGH CHILD MODELS WILL ONLY 
C-------NEED BOTTOM LAYER.  I WOULD NEED A DIFFERENT VK VARIABLE FOR 
C-------PARENT VS CHILD MODELS.
      IF(IUNIT(1) .NE. 0) THEN
        ALLOCATE(HK(NCOL,NROW,NLAY),VK(NCOL,NROW,NLAY))
      ELSE IF (IUNIT(37) .NE. 0) THEN
        ALLOCATE(HK(1,1,1),VK(NCOL,NROW,NLAY))
      ELSE IF (IUNIT(23) .NE. 0) THEN
        ALLOCATE(HK(1,1,1),VK(1,1,1))
      END IF
C
C-------SAVE POINTER DATA TO ARRARYS
      CALL SGWF2LGR2PSV(IGRID)
C
      RETURN
      END SUBROUTINE GWF2LGR2AR
      
C***********************************************************************
      
C***********************************************************************

C-----VERSION 1.0 02JULY2009 GWF2LGR2RP
      SUBROUTINE  GWF2LGR2RP(KPER,IPG,IGRID)
C     ******************************************************************
C     FIND A MAPPING BETWEEN THE COL, ROW, AND LAYER OF AN INTERFACE 
C     CHILD CELL TO THE CORRESPONDING LOCATION OF THE PARENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IBOUND,NCOL,NROW,NLAY,GLOBALDAT
      USE LGRMODULE,   ONLY:NCPP,NCPPL,NPL,NCON,KPLC,IPLC,JPLC,IEDG,
     1                      JEDG,NPCBEG,NPRBEG,NPCEND,NPREND,IBFLG,
     2                      IBOTFLG,ISCHILD,NBNODES,NPBNODES,IFACEGN,
     3                      ICBOUND,GNFLUXR
      LOGICAL DONE
C     ------------------------------------------------------------------
      CALL SGWF2LGR2PNT(IGRID)
C
C1------ONLY CREATING MAPPING IF ON THE FIRST STRESS PERIOD OF A CHILD.  
      IF(KPER .NE. 1 .OR. ISCHILD .LT. 0) RETURN
C  
C2------INITIALIZE VARIABLES
      NCON = 0
      GNFLUXR = 0.
      NK=NLAY
      IBMAX=0
      IF (IBOTFLG .EQ. 1) NK=NK-1
C3------FIND FIRST INTERFACE CELL AND STORE LOCATION
      DO I = 1,NROW
        DO J = 1,NCOL
C3A-----FIND FIRST INTERFACE CELL.  STORE LOCATION AND EXIT LOOP
          IF (IBOUND(J,I,1) .EQ. IBFLG) THEN
            JSTART=J
            ISTART=I
            GOTO  40
          ENDIF  
        ENDDO
      ENDDO
   40 CONTINUE
      JNEXT=JSTART
      INEXT=ISTART
C
C4------TRACE OUT THE INTERFACE BY MOVING IN THE CLOCKWISE DIRECTION.
C4------KEEP TRACK OF WHICH DIRECTION WE ARE MOVING AND ADJACENT 
C4------CONNECTIONS.  THE FACES AND DIRECTIONS USE MODPATH CONVENTION
C4------MOVDIR=1 LEFT; MOVDIR=2 RIGHT; MOVDIR=3 FRONT; MOVDIR=4 BACK
C4------LOOP THROUGH CHILD LAYERS.  FIND IF CURRENT CELL IS WITHIN THIS
C4------PARENT CELL.  IF SO, STORE AND EXIT

      DO K = 1, MAX(NK,1)
        KP=1
        KCCBEG=1
        KCCEND=KCCBEG+NCPPL(KP)-1
        DO KK=1,NPL
          IF(K .GE. KCCBEG .AND. K .LE. KCCEND)THEN
            KP=KK
            EXIT 
          ENDIF
          KCCBEG=KCCEND+1
          KCCEND=KCCBEG+NCPPL(KK+1)-1
        ENDDO    
        MOVDIR=4
        DONE = .FALSE.
        DO WHILE (.NOT. DONE)
        
C4A-----FIND PARENT COLUMN AND ROW LOCATIONS
          CALL SGWF2LGR2PLOC(NCPP,NPCBEG,NPCEND,J,JP)
          CALL SGWF2LGR2PLOC(NCPP,NPRBEG,NPREND,I,IP)

C4B-----FIND THE CHILD INDEX IN TERMS OF SWEEPING ALONG COL, ROW, LAYERS
          IB=0
          DO KK=1,NLAY
            DO II=1,NROW
              DO JJ=1,NCOL
                IF (ICBOUND(JJ,II,KK) .EQ. IBFLG) IB=IB+1
                IF (II .EQ. I .AND. JJ .EQ. J .AND. KK .EQ. K) GOTO 50
              ENDDO
            ENDDO
          ENDDO
  50      CONTINUE
C4C-----STORE COLUMN AND ROW OF EDGES
cswm                  write(*,*) IB, I, IP, J, JP, KP, movdir
          IF(K .EQ. 1)THEN
            JEDG(IB)=J
            IEDG(IB)=I
          ENDIF
          IBMAX=MAX(IBMAX,IB)            
C5A-D ARE SIMILAR.  
C5A-----SEARCH IN DIRECTION OF MOVEMENT.  IF NO CONNECTION IN THIS 
C5A-----DIRECTION, THEN IT IS A CORNER CELL.  DETERMINE WHETHER AN
C5A-----INSIDE OR OUTSIDE CORNER.  CHECK IF OUTSIDE OF PARENT DOMAIN
C5A-----SET GHOSTNODE FACE AND PARENT CONNECTION INDICES.
          IF(MOVDIR .EQ. 2)THEN
            IF(J+1 .LE. NCOL)THEN
              IF(ICBOUND(J+1,I,K) .EQ. IBFLG)THEN
                IF(IP-1 .GT. 0)THEN
                  NCON(IB)=1
                  IFACEGN(1,IB)=4
                  JPLC(1,IB)= JP
                  IPLC(1,IB)= IP - 1
                  KPLC(1,IB)= KP
                ENDIF
C-------FINISHED THIS CELL.  RUN THROUGH AGAIN 
                JNEXT=J+1
                GOTO 100
              ENDIF
            ENDIF
C-------CORNER CELL
            IF(ICBOUND(J,I+1,K) .EQ. IBFLG)THEN
C-------INSIDE CORNER: SET HORIZONTAL CONNECTIONS. NEW DIRECTION: FRONT
              NC=0 
              IF(IP-1 .GE. 1)THEN 
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=4
                JPLC(NC,IB)= JP
                IPLC(NC,IB)= IP - 1
                KPLC(NC,IB)= KP
              ENDIF
              IF(JP+1 .LE. GLOBALDAT(IPG)%NCOL)THEN   
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=2
                JPLC(NC,IB)= JP + 1
                IPLC(NC,IB)= IP 
                KPLC(NC,IB)= KP
              ENDIF
              INEXT=I+1
              MOVDIR=3
            ELSEIF(ICBOUND(J,I-1,K) .EQ. IBFLG)THEN
C-------OUTSIDE CORNER: NO HORIZONTAL CONNECTIONS. NEW DIRECTION: BACK
              NCON(IB)=0
              INEXT=I-1
              MOVDIR=4
            ENDIF 
C-------FINISHED THIS CELL.  RUN THROUGH AGAIN 
            GOTO 100
          ENDIF
C
C5B-----MOVE LEFT
          IF(MOVDIR .EQ. 1)THEN
            IF(J-1 .GE. 1)THEN
              IF(ICBOUND(J-1,I,K) .EQ. IBFLG)THEN
                IF(IP+1 .LE. GLOBALDAT(IPG)%NROW)THEN
                  NCON(IB)=1
                  IFACEGN(1,IB)=3
                  JPLC(1,IB)= JP
                  IPLC(1,IB)= IP + 1
                  KPLC(1,IB)= KP
                ENDIF
C-------FINISHED THIS CELL.  RUN THROUGH AGAIN 
                JNEXT=J-1
                GOTO 100
              ENDIF
            ENDIF
C-------CORNER CELL
            IF(ICBOUND(J,I-1,K) .EQ. IBFLG)THEN
C-------INSIDE CORNER: SET HORIZONTAL CONNECTIONS. NEW DIRECTION: BACK
              NC=0 
              IF(IP+1 .LE. GLOBALDAT(IPG)%NROW)THEN 
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=3
                JPLC(NC,IB)= JP
                IPLC(NC,IB)= IP + 1
                KPLC(NC,IB)= KP
              ENDIF
              IF(JP-1 .GE. 1)THEN
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=1
                JPLC(NC,IB)= JP - 1
                IPLC(NC,IB)= IP 
                KPLC(NC,IB)= KP
              ENDIF
              INEXT=I-1
              MOVDIR=4
            ELSEIF(ICBOUND(J,I+1,K) .EQ. IBFLG)THEN
C-------OUTSIDE CORNER: NO HORIZONTAL CONNECTIONS. NEW DIRECTION: FRONT
              NCON(IB)=0
              INEXT=I+1
              MOVDIR=3
            ENDIF 
C-------FINISHED THIS CELL.  SET INDICES FOR NEXT CELL RUN THROUGH AGAIN 
            GOTO 100
          ENDIF
C
C5C-----MOVE FRONTWARDS
          IF(MOVDIR .EQ. 3)THEN
            IF(I+1 .LE. NROW)THEN
              IF(ICBOUND(J,I+1,K) .EQ. IBFLG)THEN
                IF(JP+1 .LE. GLOBALDAT(IPG)%NCOL)THEN
                  NCON(IB)=1
                  IFACEGN(1,IB)=2
                  JPLC(1,IB)= JP + 1
                  IPLC(1,IB)= IP 
                  KPLC(1,IB)= KP
                ENDIF
                INEXT=I+1
C-------FINISHED THIS CELL.  RUN THROUGH AGAIN 
                GOTO 100
              ENDIF
            ENDIF
C-------CORNER CELL
            IF(ICBOUND(J-1,I,K) .EQ. IBFLG)THEN
C-------INSIDE CORNER: SET HORIZONTAL CONNECTIONS. NEW DIRECTION: LEFT
              NC=0
              IF(JP+1 .LE. GLOBALDAT(IPG)%NCOL)THEN
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=2
                JPLC(NC,IB)= JP + 1
                IPLC(NC,IB)= IP 
                KPLC(NC,IB)= KP
              ENDIF
              IF(IP+1 .LE. GLOBALDAT(IPG)%NROW)THEN
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=3
                JPLC(NC,IB)= JP 
                IPLC(NC,IB)= IP + 1
                KPLC(NC,IB)= KP
              ENDIF
              JNEXT=J-1
              MOVDIR=1
            ELSEIF(ICBOUND(J+1,I,K) .EQ. IBFLG)THEN
C-------OUTSIDE CORNER: NO HORIZONTAL CONNECTIONS. NEW DIRECTION: RIGHT
              NCON(IB)=0
              JNEXT=J+1
              MOVDIR=2
            ENDIF 
C-------FINISHED THIS CELL.  SET INDICES FOR NEXT CELL RUN THROUGH AGAIN 
            GOTO 100
          ENDIF
C
C5D-----MOVE BACKWARDS 
          IF(MOVDIR .EQ. 4)THEN
            IF(I-1 .GE. 1)THEN
              IF(ICBOUND(J,I-1,K) .EQ. IBFLG)THEN
                IF(JP-1 .GE. 1)THEN
                  NCON(IB)=1
                  IFACEGN(1,IB)=1
                  JPLC(1,IB)= JP - 1
                  IPLC(1,IB)= IP 
                  KPLC(1,IB)= KP
                ENDIF
                INEXT=I-1
C-------FINISHED THIS CELL.  RUN THROUGH AGAIN 
                GOTO 100
              ENDIF
            ENDIF
C-------CORNER CELL
            IF(ICBOUND(J+1,I,K) .EQ. IBFLG)THEN
C-------INSIDE CORNER: SET HORIZONTAL CONNECTIONS. NEW DIRECTION: RIGHT
              NC=0
              IF(JP-1 .GE. 1)THEN
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(1,IB)=1
                JPLC(NC,IB)= JP - 1
                IPLC(NC,IB)= IP 
                KPLC(NC,IB)= KP
              ENDIF
              IF(IP-1 .GE. 1)THEN
                NC=NC+1
                NCON(IB)=NC
                IFACEGN(NC,IB)=4
                JPLC(NC,IB)= JP 
                IPLC(NC,IB)= IP - 1
                KPLC(NC,IB)= KP
              ENDIF
              JNEXT=J+1
              MOVDIR=2
            ELSEIF(ICBOUND(J-1,I,K) .EQ. IBFLG)THEN
C-------OUTSIDE CORNER: NO HORIZONTAL CONNECTIONS. NEW DIRECTION: LEFT
              NCON(IB)=0
              JNEXT=J-1
              MOVDIR=1
            ENDIF 
C-------FINISHED THIS CELL.  SET INDICES FOR NEXT CELL RUN THROUGH AGAIN 
            GOTO 100
          ENDIF
  100     CONTINUE  

C-------SET INDICES FOR NEXT CELL AND CHECK IF THE TRACE IS COMPLETED 
          J=JNEXT
          I=INEXT
          IF(J .EQ. JSTART .AND. I .EQ. ISTART) DONE = .TRUE.
        END DO
C5E-----STORE THE NUMBER OF INTERFACE NODES IN THE FIRST LAYER
        IF (K.EQ.1) IB1=IBMAX
      END DO
C
            
C6------LOOP THROUGH CELLS OF THE BOTTOM LAYER.  IF IT IS AN EDGE CELL
C6------THEN ADD IN THE HORIZONTAL CONNECTIONS.  IF IT IS AN INTERIOR
C6------CELL THEN ADD BOTTOM CONNECTIONS ONLY
      IB=IBMAX
      IF(IBOTFLG .EQ. 1)THEN
        DO I=1,NROW
          DO J=1,NCOL
            IF(ICBOUND(J,I,K) .EQ. IBFLG)THEN
              IB=IB+1
              IBMAX=MAX(IBMAX,IB)
C6A-----DETERMINE IF THIS IS AN EDGE CELL.  IF SO, ADD HORIZONTAL 
C6A-----CONNECTIONS.
              DO IIB=1,IB1
                IF(JEDG(IIB) .EQ. J .AND. IEDG(IIB) .EQ. I)THEN
                  NCON(IB)=NCON(IIB)
                  DO NC=1,NCON(IIB)
                    JPLC(NC,IB)= JPLC(NC,IIB)
                    IPLC(NC,IB)= IPLC(NC,IIB) 
                    KPLC(NC,IB)= NPL 
                    IFACEGN(NC,IB) = IFACEGN(NC,IIB)
                  ENDDO
                ENDIF
              ENDDO
C6B-----FIND PARENT COLUMN AND ROW LOCATIONS OF INTERIOR CELLS
              CALL SGWF2LGR2PLOC(NCPP,NPCBEG,NPCEND,J,JP)
              CALL SGWF2LGR2PLOC(NCPP,NPRBEG,NPREND,I,IP)
C6C-----ADD CONNECTION TO BOTTOM FACE                                
              NCON(IB)=NCON(IB) + 1
              JPLC(NCON(IB),IB)= JP
              IPLC(NCON(IB),IB)= IP 
              KPLC(NCON(IB),IB)= NPL + 1
              IFACEGN(NCON(IB),IB)=6
            ENDIF
          ENDDO    
        ENDDO    
      ENDIF
      do i=1,nbnodes
        do n=1,ncon(i)
!swmdbg          write(57,*)i,n,jplc(n,i),iplc(n,i),kplc(n,i),ifacegn(n,i)
        enddo
      enddo

      IF(IBMAX .NE. NBNODES) THEN
       WRITE (*,*) 'IBMAX=',IBMAX, ' NBNODES=', NBNODES         
       CALL USTOP('ERROR IN LOCATION OF CHILD INTERFACE NODES')
      ENDIF
C        
C6------CREATE AN ORDERING FOR THE BOUNDARY NODES IN TERMS OF THE PARENT
C6------ARRAY LOOPING INDICES.  LOOP THROUGH PARENT INDICES, CHECK IF 
C6------COLUMN, ROW, AND LAYER MATCH ANY OF THE INTERFACE CELLS.  IF SO,
C6------INCREMENT INDEX COUNTER (ON FIRST OCCURENCE) AND STORE INDEX
      IBP=0 
      NPBNODES=0
      DO KP=1,GLOBALDAT(1)%NLAY
        DO IP = 1,GLOBALDAT(1)%NROW   !!!!!!! hardwired for single parent
          DO JP = 1,GLOBALDAT(1)%NCOL
!            IBSTART=0
            DO IB=1,NBNODES
              IF(JPLC(1,IB) .EQ. JP .AND. IPLC(1,IB) .EQ. IP .AND. 
     &           KPLC(1,IB) .EQ. KP) THEN 
                NPBNODES=NPBNODES+1
C--------SKIP TO NEXT JP INDEX TO AVOID REPEATS
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END SUBROUTINE  GWF2LGR2RP
      
C***********************************************************************

C-----VERSION 1.1 10OCTOBER2006 GWF2LGR1INITP
      SUBROUTINE GWF2LGR2INITP(KPER,KSTP,LGRITER,NPCBEG,NPRBEG,
     1                         NPLBEG,NPCEND,NPREND,NPLEND,ISHFLG,
     2                         MXLGRITER,INEVT,INRCH,INRES,LG,IGRID)
C     ******************************************************************
C     THIS ROUTINE ZEROS OUT THE INTERIOR OF THE PARENT GRID WHERE THE
C     THE CHILD WILL BE.  
C     NOTE: HARDWIRED FOR RECTANGULAR REFINEMENT
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NLAY,IBOUND,HNEW,CR,CC,CV
      USE GWFBASMODULE,ONLY:HNOFLO
      USE LGRMODULE,   ONLY:IBPFLG

      DOUBLE PRECISION HNF
C     ------------------------------------------------------------------
C
      ZERO = 0.
      HNF = HNOFLO
C
C1------INITIALIZE VERTICAL DISCRETIZATION LOOP LIMITS
      IF(NPLBEG .EQ. NPLEND)THEN
        KBEG=NPLBEG
        KEND=KBEG
      ELSE
        KBEG=1
        KEND=NPLEND
      ENDIF
C
C2------CALLED AFTER THE FIRST ITERATION IF NOT RUNNING 1-WAY COUPLED
C2------AFTER THE FIRST ITERATION OF THE FIRST TIME STEP OF THE FIRST 
C2------STRESS  PERIOD, ZERO OUT THE INTERIOR.
C2------REWETTING MIGHT CAUSE CELLS TO REACTIVATE, SO NEED TO RE-ZERO.
      IF(MXLGRITER .GT. 1 .AND. IGRID .EQ. 1) THEN
        DO K = KBEG,KEND
! should check laywet and kper kstep lgriter to see if this only needs
        ! to be done once.
          DO I = NPRBEG,NPREND
            DO J= NPCBEG,NPCEND
cswm              IF(I .EQ. NPRBEG .OR. I .EQ. NPREND .OR. 
cswm     &           J .EQ. NPCBEG .OR. J .EQ. NPCEND)THEN
cswm                IBOUND(J,I,K) = IBPFLG(LG) 
C2A-----PICK UP BOTTOM  
cswm              ELSEIF(K .EQ. KEND .AND. IBOTFLG .EQ. 1)THEN
cswm                IBOUND(J,I,K) = IBPFLG(LG)
cswm              ELSE
C2B-----ZERO OUT INTERIOR
                IBOUND(J,I,K)=0
                IF(ISHFLG .EQ. 0 .OR. ISHFLG .EQ.1 .AND. LGRITER .GT.1) 
     1            HNEW(J,I,K) = HNF
                IF(K.NE.NLAY) CV(J,I,K)=ZERO
                IF(K.NE.1) CV(J,I,K-1)=ZERO
                CC(J,I,K)=ZERO
                IF(I.NE.1) CC(J,I-1,K)=ZERO
                CR(J,I,K)=ZERO
                IF(J.NE.1) CR(J-1,I,K)=ZERO
cswm              ENDIF
C2C-----ZERO OUT STRESS PACKAGE VALUES THAT HAVE OPTION 3
              IF(INEVT .NE. 0) CALL SGWF2LGR2EVT(J,I)
              IF(INRCH .NE. 0) CALL SGWF2LGR2RCH(J,I)
              IF(INRES .NE. 0) CALL SGWF2LGR2RES(J,I)
            ENDDO
          ENDDO
        ENDDO
C3
C3------SET HNEW TO HNOFLO ON SECOND LGR ITERATION. NOT OVERWRITTEN ON
C3------FIRST ITERATION BECAUSE HNEW WAS USED FOR CHILD START HEAD.
cswm      ELSEIF(ISHFLG .EQ. 1 .AND. LGRITER .EQ. 2 .AND. IGRID .EQ. 1 .AND.
cswm     &       KPER .EQ. 1 .AND. KSTP .EQ. 1)THEN
cswm       DO K = KBEG,KEND
cswm          DO I = NPRBEG,NPREND
cswm            DO J= NPCBEG,NPCEND
cswm              IF(IBOUND(J,I,K) .EQ. 0) HNEW(J,I,K) = HNF 
cswm            ENDDO
cswm          ENDDO
cswm        ENDDO
      ENDIF     
C4
      RETURN
      END SUBROUTINE GWF2LGR2INITP
      
C***********************************************************************

C-----VERSION 1.2 06JANUARY2009 GWF2LGR2PFM      
      SUBROUTINE GWF2LGR2PFM(KPER,KSTP,KITER,LGRITER,IUBCF,IULPF,IUHUF,
     1                       LG)
C
C     ******************************************************************
C     ADD GHOST-NODE CONDUCTANCE TERMS (CALCULATED FROM PREVIOUS CHILD?)  
C     TO HCOF AND RHS OF PARENT GRID   (CURRENTLY UPDATING W/ NEW GNCOND
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IBOUND,HNEW,BOTM,DELR,DELC,HCOF,RHS,LAYHDT,
     1                      GLOBALDAT,nlay,nrow,ncol
      USE LGRMODULE,   ONLY:IBPFLG,LGRDAT
C     ------------------------------------------------------------------
C
C1------IF ON THE VERY FIRST ITERATION, RETURN
C1------THE FIRST LGR ITERATION OF THE FIRST TIME STEP OF THE FIRST
C1------STRESS PERIOD IS A FULL PARENT GRID
      IF(LGRITER .EQ. 1 .AND. KSTP .EQ. 1 .AND. KPER .EQ. 1) RETURN
C2------ADJUST HCOF AND RHS OF PARENT FOR THE GHOST-NODE CONNECTIONS 
C
      IUBCFC=GLOBALDAT(LG)%IUNIT(1)
      IULPFC=GLOBALDAT(LG)%IUNIT(23)
      IUHUFC=GLOBALDAT(LG)%IUNIT(37)
      
C2-----LOOP THROUGH CHILD LAYERS, ROWS, AND COLUMNS AND FIND INTERFACE
C2-----SKIP IF INACTIVE
      IB=0
      DO K=1,GLOBALDAT(LG)%NLAY
        LAYHDTC=GLOBALDAT(LG)%LAYHDT(K)
        DO I=1,GLOBALDAT(LG)%NROW
          DO J=1,GLOBALDAT(LG)%NCOL
            IF (LGRDAT(LG)%ICBOUND(J,I,K) .EQ. IBPFLG(LG)) THEN
              IB=IB+1            
              IF (GLOBALDAT(LG)%IBOUND(J,I,K) .EQ. 0) THEN
                LGRDAT(LG)%GNCOND(:,IB) = 0. 
                CYCLE
              END IF
              
C2A-----FOR EACH CHILD INTERFACE CELL CALCULATE CELL THICKNESS
              IF (LAYHDTC .EQ. 0) THEN
                THICKC=GLOBALDAT(LG)%BOTM(J,I,K-1) - 
     1                 GLOBALDAT(LG)%BOTM(J,I,K)
              ELSE
                THICKC=GLOBALDAT(LG)%HNEW(J,I,K) - 
     1                 GLOBALDAT(LG)%BOTM(J,I,K)  
              END IF

C3------FOR EACH GHOST-NODE CONNECTION, CALCULATE THE GNCOND
C3------AND ADJUST HCOF AND RHS.  SKIP IF CELL IS INACTIVE
              DO NC=1,LGRDAT(LG)%NCON(IB)
C3A-----SET PARENT CELL INDICES. SKIP IF PARENT IS INACTIVE
                KP=LGRDAT(LG)%KPLC(NC,IB)
                IP=LGRDAT(LG)%IPLC(NC,IB)
                JP=LGRDAT(LG)%JPLC(NC,IB)
                IF (IBOUND(JP,IP,KP) .EQ. 0) THEN
                  LGRDAT(LG)%GNCOND(NC,IB) = 0. 
                  CYCLE
                END IF
C3B-----CALCULATE PARENT CELL DIMENSIONS
                IF (LAYHDT(KP) .EQ. 0) THEN
                  THICKP=BOTM(JP,IP,KP-1) - BOTM(JP,IP,KP)
                  LAYHDTFLG=0
                ELSE
                  THICKP=HNEW(JP,IP,KP) - BOTM(JP,IP,KP)
                  THICKPC=HNEW(JP,IP,KP) - GLOBALDAT(LG)%BOTM(J,I,K)
                  LAYHDTFLG=1
                END IF
                DELRP = DELR(JP)
                DELCP = DELC(IP)
C
C3C-----CALCULATE AND STORE GHOST-NODE CONDUCTANCE (CGN) IF FIRST LGR
C3C-----ITERATION OR SATURATED THICKNESS IS HEAD DEPENDENT 
                IF (KITER .EQ. 1 .AND. LGRITER .EQ. 1 .OR. 
     1              LAYHDTC .NE. 0 .OR. LAYHDT(KP) .NE. 0) THEN 
cswm                  LGRDAT(LG)%GNCOND=0.
                  CALL SGWF2LGR2GNCOND(I,J,K,JP,IP,KP,IUBCF,IULPF,
     1                                  IUHUF,IUBCFC,IULPFC,IUHUFC,
     2                                  IB,NC,LAYHDTFLG,THICKP,THICKPC,
     3                                  DELRP,DELCP,THICKC,LG,CGN)
C
                  LGRDAT(LG)%GNCOND(NC,IB) = CGN
                END IF
C
C3E-----ADJUST HCOF AND RHS 
                IF(LGRDAT(LG)%RELAXF .EQ. 1.) THEN
                  HCOF(JP,IP,KP)=HCOF(JP,IP,KP)-LGRDAT(LG)%GNCOND(NC,IB)
                  RHS(JP,IP,KP)=RHS(JP,IP,KP)+LGRDAT(LG)%GNCOND(NC,IB)* 
     1                        (LGRDAT(LG)%DHGN(NC,IB) -
     2                         GLOBALDAT(LG)%HNEW(J,I,K))
                ELSE
                  RHS(JP,IP,KP)=RHS(JP,IP,KP)+LGRDAT(LG)%GNFLUXR(NC,IB) 
                ENDIF
                   

!swmdbg       If(kiter .eq. 1 .and. lgriter .eq. 2)   
!swmdbg     &                    write(51,*) ib,k,i,j,nc,kp,ip,jp,
!swmdbg     &                     LGRDAT(LG)%IFACEGN(NC,IB),
!swmdbg     1                     LGRDAT(LG)%GNCOND(NC,IB),
!swmdbg     2                     LGRDAT(LG)%DHGN(NC,IB),
!swmdbg     3                     GLOBALDAT(LG)%HNEW(J,I,K),HCOF(JP,IP,KP),
!swmdbg     4                     RHS(JP,IP,KP)


              ENDDO     
            END IF
          ENDDO     
        ENDDO     
      ENDDO     
C4
      do k=1,nlay
!swmdbg        write(61,*) 'nlay=',k
!swmdbg        write(62,*) 'nlay=',k
      do i=1,nrow
!swmdbg         write(61,'(100f11.5)') (hcof(j,i,k),j=1,ncol)
!swmdbg         write(62,'(100f11.5)') (rhs(j,i,k),j=1,ncol)
         enddo
         enddo
      RETURN
      END SUBROUTINE GWF2LGR2PFM
      
C***********************************************************************

C-----VERSION 1.2 06JANUARY2009 GWF2LGR2CFM      
      SUBROUTINE GWF2LGR2CFM(KITER,LGRITER,IUBCF,IULPF,IUHUF,
     1                       LG)
C
C     ******************************************************************
C     ADD GHOST-NODE CONDUCTANCE TERMS TO HCOF AND RHS OF CHILD GRID
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,BOTM,HCOF,RHS,
     &                      LAYHDT,GLOBALDAT
      USE LGRMODULE,   ONLY:NCON,JPLC,IPLC,KPLC,IBFLG,ICBOUND,GNHEAD,
     &                      GNCOND
C
C     ------------------------------------------------------------------
C
C1------STORE UNIT NUMBER FOR PARENT FLOW PACKAGE
Cswm: NOTE: hardwired for single parent grid 
      IUBCFP=GLOBALDAT(1)%IUNIT(1)
      IULPFP=GLOBALDAT(1)%IUNIT(23)
      IUHUFP=GLOBALDAT(1)%IUNIT(37)

C
C3-----LOOP THROUGH CHILD LAYERS, ROWS, AND COLUMNS AND FIND INTERFACE
C3-----SKIP IF INACTIVE
      IB=0
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
            IF (ICBOUND(J,I,K) .EQ. IBFLG) THEN
              IB=IB+1            
              IF (IBOUND(J,I,K) .EQ. 0) THEN
                GNCOND(:,IB) = 0.
                CYCLE
              END IF
C3A-----FOR EACH CHILD INTERFACE CELL CALCULATE CELL THICKNESS
              IF (LAYHDT(K) .EQ. 0) THEN
                THICKC=BOTM(J,I,K-1) - BOTM(J,I,K)
              ELSE
                THICKC=HNEW(J,I,K) - BOTM(J,I,K)  
              END IF
C
C4------FOR EACH GHOST-NODE CONNECTION, CALCULATE THE GNCOND
C4------AND ADJUST HCOF AND RHS. SKIP IF INACTIVE
              DO NC=1,NCON(IB)
C4A-----SET PARENT CELL INDICES. SKIP IF PARENT IS INACTIVE
                KP=KPLC(NC,IB)
                IP=IPLC(NC,IB)
                JP=JPLC(NC,IB)
                IF (GLOBALDAT(1)%IBOUND(JP,IP,KP) .EQ. 0) THEN
                  GNCOND(NC,IB) = 0.
                  CYCLE
                END IF
C4B-----CALCULATE PARENT CELL DIMENSIONS
Cswm: NOTE: hardwired for single parent grid
                IF (GLOBALDAT(1)%LAYHDT(KP) .EQ. 0) THEN
                  THICKP=GLOBALDAT(1)%BOTM(JP,IP,KP-1) - 
     1                   GLOBALDAT(1)%BOTM(JP,IP,KP)
                  LAYHDTFLG = 0
                ELSE
                  THICKP=GLOBALDAT(1)%HNEW(JP,IP,KP) - 
     1                   GLOBALDAT(1)%BOTM(JP,IP,KP)
                  THICKPC=GLOBALDAT(1)%HNEW(JP,IP,KP) - BOTM(J,I,K)
                  LAYHDTFLG = 1
                END IF
                DELRP = GLOBALDAT(1)%DELR(JP)
                DELCP = GLOBALDAT(1)%DELC(IP)
      
C4C-----CALCULATE AND STORE GHOST-NODE CONDUCTANCE (CGN) IF FIRST LGR
C4C-----ITERATION OR SATURATED THICKNESS IS HEAD DEPENDENT 
                IF (KITER .EQ. 1 .AND. LGRITER .EQ. 1 .OR. 
     1              LAYHDT(K) .NE. 0 .OR. 
     2              GLOBALDAT(1)%LAYHDT(KP) .NE. 0) THEN 
                  CALL SGWF2LGR2GNCOND(I,J,K,JP,IP,KP,IUBCFP,IULPFP,
     1                                  IUHUFP,IUBCF,IULPF,IUHUF,
     2                                  IB,NC,LAYHDTFLG,THICKP,THICKPC,
     3                                  DELRP,DELCP,THICKC,LG,CGN)
                  GNCOND(NC,IB) = CGN
                END IF
C
C4D-----ADJUST HCOF AND RHS 
                  HCOF(J,I,K)=HCOF(J,I,K) - GNCOND(NC,IB)
                  RHS(J,I,K)=RHS(J,I,K) - GNCOND(NC,IB)*GNHEAD(NC,IB)

              ENDDO     
            END IF
          ENDDO     
        ENDDO     
      ENDDO     
C5------RETURN
      RETURN
      END SUBROUTINE GWF2LGR2CFM
      
c***********************************************************************

*      SUBROUTINE GWF2LGR2DARCY(KPER,KSTP,LGRITER,IPG,LG)
*C     ******************************************************************
*C     THIS ROUTINE USES A DARCY-PLANAR INTERPOLATION BETWEEN THE PARENT 
*C     HEADS ONTO THE GHOST NODES. 
*C     IF ISHFLG .NE. 0 AND ON THE FIRST ITERATION OF THE FIRST TIME STEP 
*C     OF THE FIRST STRESS PERIOD, THEN ALSO USE THE PARENT SOLUTION TO 
*C     INTIALIZE THE INTERIOR SOLUTION OF THE CHILD GRID.
*C
*C     NOTE: NEED TO ADD CHECKS FOR REFINEMENT BEGINNING/ENDING AT EDGE
*C     OF PARENT GRID (DOES IT MATTER?  MIGHT BE EASIER TO SET GNCOND=0)
*C
*C        SPECIFICATIONS:
*C     ------------------------------------------------------------------
*      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,HNEW,DELR,DELC,BOTM,HOLD,
*     1                      STRT,GLOBALDAT
*      USE LGRMODULE,   ONLY:NPCBEG,NPRBEG,NPCEND,NPREND,NPLEND,NCPP,NPL,
*     1                      NCPPL,ISHFLG,IBOTFLG,IBFLG,MXLGRITER,IBMAXH,
*     2                      NCMAXH,RELAXH,HDIFFM,NODEH,NCON,JPLC,IPLC,
*     3                      KPLC,IFACEGN,ICBOUND,GNHEAD,DHGN 
*      DOUBLE PRECISION      HP, HP2
*C     ------------------------------------------------------------------
*C
*C
*C1------INITIALIZE VARIABLES OF THE PARENT GRID
*      IUBCFP=GLOBALDAT(IPG)%IUNIT(1)
*      IULPFP=GLOBALDAT(IPG)%IUNIT(23)
*      IUHUFP=GLOBALDAT(IPG)%IUNIT(37)
*      NROWP=GLOBALDAT(IPG)%NROW
*      NCOLP=GLOBALDAT(IPG)%NCOL
*      NLAYP=GLOBALDAT(IPG)%NLAY
*      ICEND=NCPP
*      JCEND=NCPP
*
*C2------INITIALIZE MAXIMUM HEAD DIFFERENCE AND ITS LOCATION
*      HDIFFM=0.
*      NODEH(1)=1
*      NODEH(2)=1
*      NODEH(3)=1
*      IBMAXH=1
*      NCMAXH=1
*      
*      !swm: make this a subroutine.
*C2------IF ON THE VERY FIRST ITERATION AND IF ISHFLG IS SET TO 1
*C2------INITIALIZE WITH PARENT HEADS 
*      IF(ISHFLG .EQ. 1 .AND. KPER .EQ. 1 .AND. KSTP .EQ. 1 .AND. 
*     &   LGRITER .EQ. 1)THEN
*C2A-----BIN PARENT HEADS ONTO THE ENTIRE CHILD GRID
*        KCLAY=1
*        DO KPLAY = 1,NPL
*          KCEND=NCPPL(KPLAY)       
*          ICR=1
*          DO IPR = NPRBEG,NPREND
*            JCC=1
*            DO JPC = NPCBEG,NPCEND
*              DO K=1,KCEND
*                KM1=K-1
*                DO I=1,ICEND
*                  IM1=I-1
*                  DO J=1,JCEND
*                    JM1=J-1
*C-------BE CAREFUL NOT TO OVERWRITE CONSTANT HEADS UNLESS THEY ARE 
*C-------ALONG THE BOUNDARY.  ALSO OVERWRITE STRT AND HOLD WITH HNEW 
*C-------SO CONSISTENT VALUES ARE USED FOR INITIAL CONDITIONS WHEN 
*C-------ISHFLG=1 AND SIMULATIONS ARE TRANSIENT.
*                    IF(ICBOUND(JCC+JM1,ICR+IM1,KCLAY+KM1) .GT. 0 .OR.
*     &                 ABS(ICBOUND(JCC+JM1,ICR+IM1,KCLAY+KM1)).EQ.IBFLG)
*     &              HNEW(JCC+JM1,ICR+IM1,KCLAY+KM1)=
*     &              GLOBALDAT(IPG)%HNEW(JPC,IPR,KPLAY)
*                    STRT(JCC+JM1,ICR+IM1,KCLAY+KM1)=
*     &              HNEW(JCC+JM1,ICR+IM1,KCLAY+KM1)
*                    HOLD(JCC+JM1,ICR+IM1,KCLAY+KM1)=
*     &              HNEW(JCC+JM1,ICR+IM1,KCLAY+KM1)     
*                  ENDDO
*                ENDDO
*              ENDDO
*
*              JCC=JCC+JCEND
*            ENDDO 
*            ICR=ICR+ICEND
*          ENDDO
*          KCLAY=KCLAY+KCEND
*        ENDDO
*C
*C3------NOT ON FIRST ITERATION, OR RUNNING 1-WAY COUPLED, SO ONLY 
*C3------IF 1-WAY COUPLED OR FIRST LGR ITERATION, NO RELAXATION
*      END IF
*      TMP=RELAXH
*      IF(MXLGRITER .LE. 1 .OR. (LGRITER .EQ. 1 .AND. KSTP .EQ. 1
*     &                          .AND. KPER .EQ. 1)) THEN
*         RELAXH=1.0
*         DHGN = 0.
*      ENDIF
*C   
*      IB=0
*      JP0 = JPLC(1,1)
*      IP0 = IPLC(1,1)
*      KP0 = KPLC(1,1)
*      IF (IFACEGN(1,1) .EQ. 3) JP0=JP0-1
*      IF (IFACEGN(1,1) .EQ. 1) IP0=IP0-1
*      DPX=0.
*      DPY=0.
*      XP0=0.
*      YP0=0.
*      IF(JP0 .GE. 1 .AND. IP0 .GE. 1) THEN
*        CALL SGWF2LGR2XYZLOC(JP0,IP0,KP0,NCOLP,NROWP,NLAYP,
*     1                      GLOBALDAT(IPG)%DELR,GLOBALDAT(IPG)%DELC,
*     2                      GLOBALDAT(IPG)%BOTM,XP0,YP0,ZP0)
*        DPX = 0.5*GLOBALDAT(IPG)%DELR(JP0)
*        DPY = 0.5*GLOBALDAT(IPG)%DELC(IP0)
*      ELSEIF(JP0 .GE. 1 .AND. IP0 .LT. 1) THEN 
*        CALL SGWF2LGR2XYZLOC(JP0,1,KP0,NCOLP,NROWP,NLAYP,
*     1                      GLOBALDAT(IPG)%DELR,GLOBALDAT(IPG)%DELC,
*     2                      GLOBALDAT(IPG)%BOTM,XP0,YP0,ZP0)
*        DPX = 0.5*GLOBALDAT(IPG)%DELR(JP0)
*        YP0 = 0.
*      ELSEIF(JP0 .LT. 1 .AND. IP0 .GE. 1) THEN 
*        CALL SGWF2LGR2XYZLOC(1,IP0,KP0,NCOLP,NROWP,NLAYP,
*     1                      GLOBALDAT(IPG)%DELR,GLOBALDAT(IPG)%DELC,
*     2                      GLOBALDAT(IPG)%BOTM,XP0,YP0,ZP0)
*        DPY = 0.5*GLOBALDAT(IPG)%DELC(IP0)
*        XP0 = 0.
*      ENDIF
*
*        do k=1,nlayp
*!swmdbg          write(53,*) 'nlayp=',k
*        do i=1,globaldat(ipg)%nrow
*!swmdbg      write(53,'(100f11.5)') (globaldat(ipg)%hnew(j,i,k),
*!swmdbg     1                        j=1,globaldat(ipg)%ncol)
*        enddo
*        enddo
*      DO K = 1,NLAY
*        DO I = 1,NROW
*          DO J = 1,NCOL 
*            IF(ICBOUND(J,I,K) .EQ. IBFLG)THEN
*              IB=IB+1
*              CALL SGWF2LGR2XYZLOC(J,I,K,NCOL,NROW,NLAY,DELR,DELC,BOTM,
*     1                            XC1,YC1,ZC1)
*              DO NC = 1, NCON(IB)
*                JP= JPLC(NC,IB)
*                IP= IPLC(NC,IB)
*                KP= KPLC(NC,IB)
*C-------IF PARENT CELL IS INACTIVE, THEN SKIP
*                IF(GLOBALDAT(IPG)%IBOUND(JP,IP,KP) .EQ. 0) CYCLE
*                HP = GLOBALDAT(IPG)%HNEW(JP,IP,KP)
*                IFACE =IFACEGN(NC,IB)
*                IF (GLOBALDAT(IPG)%LAYHDT(KP) .EQ. 0) THEN
*                  THICKP=GLOBALDAT(IPG)%BOTM(JP,IP,KP-1) -
*     1                   GLOBALDAT(IPG)%BOTM(JP,IP,KP)
*                ELSE
*                  THICKP=GLOBALDAT(IPG)%HNEW(JP,IP,KP) -
*     1                   GLOBALDAT(IPG)%BOTM(JP,IP,KP)
*                ENDIF
*
*                CALL SGWF2LGR2HKP(JP,IP,KP,IUBCFP,IULPFP,IUHUFP,HKPX,
*     1                            HKPY,HKPZ,IPG)
*                CALL SGWF2LGR2XYZLOC(JP,IP,KP,NCOLP,NROWP,NLAYP,
*     1                              GLOBALDAT(IPG)%DELR,
*     2                              GLOBALDAT(IPG)%DELC,
*     3                              GLOBALDAT(IPG)%BOTM,XP,YP,ZP)
*                DHX=0.
*                DHY=0.
*                DHZ=0.
*
*C-------GET X COMPONENT
*                IF (IFACE .NE. 1 .AND. IFACE .NE. 2) THEN
*                  XC = XC1 + DPX
*                  XP = XP - XP0
*                  DX=XP-XC
*                  AREA = THICKP*GLOBALDAT(IPG)%DELC(IP)
*                  
*                  IF(DX .GT. 0. .AND. JP-1 .GE. 1) THEN
*                      
*                    HP2=GLOBALDAT(IPG)%HNEW(JP-1,IP,KP)
*                    CRP = GLOBALDAT(IPG)%CR(JP-1,IP,KP) 
*                    IF (GLOBALDAT(IPG)%IBOUND(JP-1,IP,KP) .EQ. 0) 
*     1               CRP = 0.
*!swmdbg                    DXP= 0.5*(GLOBALDAT(IPG)%DELR(JP) +  
*!swmdbg     1                        GLOBALDAT(IPG)%DELR(JP-1))
*         jp2=jp-1
*                    DHX = -CRP/HKPX/AREA*DX*(HP-HP2)
*                  ELSEIF(DX .LT. 0. .AND. JP+1 .LE. NCOLP) THEN  
*                    HP2=GLOBALDAT(IPG)%HNEW(JP+1,IP,KP)
*                    CRP = GLOBALDAT(IPG)%CR(JP,IP,KP) 
*                    IF (GLOBALDAT(IPG)%IBOUND(JP+1,IP,KP) .EQ. 0) 
*     1               CRP = 0.
*!swmdbg                    DXP=  0.5*(GLOBALDAT(IPG)%DELR(JP) + 
*!swmdbg     1                         GLOBALDAT(IPG)%DELR(JP+1))
*                    DX=-DX
*         jp2=jp+1
*                    DHX = -CRP/HKPX/AREA*DX*(HP-HP2)
*                  ENDIF
*
*!swmdbg      write(53,*)ib,k,i,j,nc,iface,kp,ip,jp,xc,xp,DX,DXP,HP,HP2,DHX,
*!swmdbg     1           jp2,CRP,HKPX,AREA
*                END IF
*C-------GET Y COMPONENT
*                IF (IFACE .NE. 3 .AND. IFACE .NE. 4) THEN
*                  YC = YC1 + DPY
*                  YP = YP - YP0
*                  DY=YP-YC
*                  AREA = THICKP*GLOBALDAT(IPG)%DELR(JP)
*                  IF(DY .GT. 0. .AND. IP-1 .GE. 1) THEN
*                    HP2=GLOBALDAT(IPG)%HNEW(JP,IP-1,KP)
*                    CCP = GLOBALDAT(IPG)%CC(JP,IP-1,KP) 
*                    IF (GLOBALDAT(IPG)%IBOUND(JP,IP-1,KP) .EQ. 0) 
*     1               CCP = 0.
*!swmdbg                    DYP= 0.5*(GLOBALDAT(IPG)%DELC(IP) +
*!swmdbg     1                        GLOBALDAT(IPG)%DELC(IP-1))
*         ip2=ip-1
*
*                    DHY = -CCP/HKPY/AREA*DY*(HP-HP2)
*                  ELSEIF(DY .LT. 0. .AND. IP+1 .LE. NROWP) THEN
*                    HP2=GLOBALDAT(IPG)%HNEW(JP,IP+1,KP)
*                    CCP = GLOBALDAT(IPG)%CC(JP,IP,KP) 
*                    IF (GLOBALDAT(IPG)%IBOUND(JP,IP+1,KP) .EQ. 0) 
*     1               CCP = 0.
*!swmdbg                    DYP= 0.5*(GLOBALDAT(IPG)%DELC(IP) +
*!swmdbg     1                        GLOBALDAT(IPG)%DELC(IP+1))
*                    DY=-DY
*         ip2=ip+1
*                    DHY = -CCP/HKPY/AREA*DY*(HP-HP2)
*                  ENDIF
*!swmdbg      write(53,*)ib,k,i,j,nc,iface,kp,ip,jp,yc,yp,DY,DYP,HP,HP2,DHY,
*!swmdbg     1             ip2,CCP,HKPY,AREA
*                END IF
*                IF (IFACE .NE. 5 .AND. IFACE .NE. 6) THEN
*C-------GET Z COMPONENT
*                  DZ=ZP-ZC1
*cswm                  IF(NPL .GT. 1 .AND. (DZ .GT. 0. .OR. KP .GT. 1)) THEN
*                    AREA=GLOBALDAT(IPG)%DELC(IP)*GLOBALDAT(IPG)%DELR(JP)
*
*                    IF(DZ .GT. 0. .AND. KP+1 .LE. NLAYP) THEN
*                      HP2=GLOBALDAT(IPG)%HNEW(JP,IP,KP+1)
*                      CVP = GLOBALDAT(IPG)%CV(JP,IP,KP) 
*                      IF (GLOBALDAT(IPG)%IBOUND(JP,IP,KP+1) .EQ. 0) 
*     1                 CVP = 0.
*!swmdbg                      DZP=  0.5*(GLOBALDAT(IPG)%BOTM(JP,IP,KP-1) - 
*!swmdbg     1                         GLOBALDAT(IPG)%BOTM(JP,IP,KP+1))
*                      DHZ = -CVP/HKPZ/AREA*DZ*(HP-HP2)
*                    ELSEIF(DZ .LT. 0. .AND. KP-1 .GE. 1) THEN
*                      HP2=GLOBALDAT(IPG)%HNEW(JP,IP,KP-1)
*                      CVP = GLOBALDAT(IPG)%CV(JP,IP,KP-1) 
*                      IF (GLOBALDAT(IPG)%IBOUND(JP,IP,KP-1) .EQ. 0) 
*     1                 CVP = 0.
*!swmdbg                      DZP= 0.5*(GLOBALDAT(IPG)%BOTM(JP,IP,KP-2) -  
*!swmdbg     1                        GLOBALDAT(IPG)%BOTM(JP,IP,KP))
*                      DZ=-DZ
*                      DHZ = -CVP/HKPZ/AREA*DZ*(HP-HP2)
*                    ENDIF
*!swmdbg      write(53,*)ib,k,i,j,nc,iface,kp,ip,jp,zc1,zp,DZ,DZP,HP,HP2,DHZ,
*!swmdbg     1           jp2,CVP,HKPZ,AREA
*cswm                  ENDIF
*                ENDIF
*C-------ADD UP X,Y, AND Z COMPONENTS TO GET CHANGE IN GHOST-NODE HEAD
*                DHGN(NC,IB) = DHX + DHY + DHZ
*
*C------FIND THE LOCATION OF THE MAXIMUM GHOST-NODE HEAD CHANGE AND STORE
*                GNHEADOLD=GNHEAD(NC,IB)
*                GNHEAD(NC,IB) = HP + DHGN(NC,IB)
*                GNHEAD(NC,IB) = RELAXH*GNHEAD(NC,IB) + 
*     1                          (1.0-RELAXH)*GNHEADOLD
*                HDIFF = GNHEAD(NC,IB) - GNHEADOLD
*                IF(ABS(HDIFF) .GT. ABS(HDIFFM))THEN
*                  HDIFFM=HDIFF
*                  NODEH(1)=K
*                  NODEH(2)=I
*                  NODEH(3)=J
*                  IBMAXH=IB
*                  NCMAXH=NC
*                ENDIF
*                        
*              END DO
*            END IF
*           END DO
*         END DO
*       END DO
*C
*C4 
*      RELAXH=TMP
*C5
*      RETURN
*      END SUBROUTINE GWF2LGR2DARCY
      
C***********************************************************************

C-----VERSION 1.0 18JULY2008 GWF2LGR2FMBF     
      SUBROUTINE GWF2LGR2FMBF(KSTP,KPER,LGRITER)
C     ******************************************************************
C     CALCULATE THE CHILD INTERFACE GHOST-NODE FLUXES 
C     POSITIVE IS INFLOW TO CHILD; NEGATIVE IS OUTFLOW FROM THE CHILD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,GLOBALDAT
      USE LGRMODULE,   ONLY:IBFLG,MXLGRITER,IBMAXF,NCMAXF,RELAXF,FDIFFM,
     1                      NODEF,JPLC,IPLC,KPLC,NCON,GNHEAD,ICBOUND,
     2                      GNFLUX,GNFLUXR,GNFLUXOLD,GNCOND  
C     ------------------------------------------------------------------

C1------IF RUNNING 1-WAY COUPLED, DO NOT CALCULATE FLUXES
cswm: keep these to get proper budget for 1way coupled      IF(MXLGRITER .LE. 1) RETURN   
C
C2------INITIALIZE GHOSTNODE FLUXES
C2------INITIALIZE MAXIMUM FLUX DIFFERENCE AND ITS LOCATION,  AND SET 
C2------RELAXATION FOR FIRST ITERATION
      GNFLUXOLD = 0.
      GNFLUX = 0.
      FDIFFM=0.
      NODEF(1)=1
      NODEF(2)=1
      NODEF(3)=1
      IBMAXF=1
      NCMAXF=1
      TMP=RELAXF
      IF(LGRITER .EQ. 1 .AND. KSTP .EQ. 1 .AND. KPER .EQ.1) RELAXF=1.0

C3------LOOP THROUGH EACH CELL AND CALCULATE FLOW INTO MODEL FROM EACH
C3------INTERFACE CELL (IBOUND .EQ. IBFLG)
      IB = 0
      DO K=1,NLAY
        DO I=1,NROW
          DO J=1,NCOL
C
C4------IF CELL IS NOT AN INTERFACE CELL SKIP IT & GO ON TO NEXT CELL.
            IF (ICBOUND(J,I,K) .NE. IBFLG) CYCLE
            IB = IB + 1
C4A-----IF AN INTERFACE CELL BUT INACTIVE, SKIP IT & GO ON TO NEXT CELL.
            IF (IBOUND(J,I,K) .EQ. 0) CYCLE
C
C5------FOR EACH INTERFACE CELL, LOOP THROUGH ALL GHOST-NODE CONNECTIONS
C5------SKIP IF INACTIVE
            DO NC=1,NCON(IB)
              JP= JPLC(NC,IB)
              IP= IPLC(NC,IB)
              KP= KPLC(NC,IB)
              IF (GLOBALDAT(1)%IBOUND(JP,IP,KP) .EQ. 0) CYCLE
C6------CALCULATE FLUX THROUGH THE GHOST-NODE. APPLY RELAXATION.
C6------STORE LOCATION OF MAXIMUM FLUX CHANGE. 
C6------GNFLUXOLD IS PREV RELAXED, GNFLUX IS CURRENT, GNFLUXR IS RELAXED 
              GNFLUXOLD(NC,IB)=GNFLUXR(NC,IB)
              GNFLUX(NC,IB) = GNCOND(NC,IB)*
     1         (GNHEAD(NC,IB) - HNEW(J,I,K))
              GNFLUXR(NC,IB) = RELAXF*GNFLUX(NC,IB) +
     1         (1.0-RELAXF)*GNFLUXOLD(NC,IB)
C
              FDIFF = (GNFLUX(NC,IB) - GNFLUXOLD(NC,IB))/
     1                 MAX(ABS(GNFLUXOLD(NC,IB)),1.0)
              IF(ABS(FDIFF) .GT. ABS(FDIFFM))THEN
                FDIFFM=FDIFF
                NODEF(1)=K
                NODEF(2)=I
                NODEF(3)=J
                IBMAXF=IB
                NCMAXF=NC
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C7------RESTORE ORIGINAL RELAXATION
      RELAXF=TMP
C8
      RETURN
      END SUBROUTINE GWF2LGR2FMBF
      
C***********************************************************************

C-----VERSION 1.1 13OCTOBER2006 GWF2LGR2CNVG
      SUBROUTINE GWF2LGR2CNVG(IGRID,NGRIDS,LGRCNVG,LGRITER,KPER,KSTP)
C     ******************************************************************
C     CHECK CONVERGENCE OF LGR
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT
      USE LGRMODULE,   ONLY:MXLGRITER,IOUTLGR,HCLOSELGR,FCLOSELGR,
     1                      GNFLUX,GNHEAD,HDIFFM,FDIFFM,NODEH,NODEF,
     2                      IBMAXH,IBMAXF,NCMAXH,NCMAXF
      use GwfBasModule, only: SGWF2BAS7PNT
C     ------------------------------------------------------------------
  500 FORMAT(1X,'LGRITER= ',I3)
  501 FORMAT(1X,'GRID NUMBER= ',I3)
  505 FORMAT(1X,'HCLOSELGR= ',1P,E10.3,' HDIFFM= ', E10.3,2X,
     &          '(',I0,',',I0,',',I0,')',' GNHEAD= ', E10.3)
  510 FORMAT(1X,'FCLOSELGR= ',1PE10.3,' FDIFFM= ', E10.3,2X,
     &          '(',I0,',',I0,',',I0,')',' GNFLUX= ',E10.3,/)
  515 FORMAT(1X,'MXLGRITER EXCEEDED FOR GRID NUMBER',I3,/,
     &          'CHECK LGR BUDGET OF G-N FLUXES TO ASSESS QUALITY OF ',
     &          'THE LGR SOLUTION',/)
C1------CHECK IF ALL GRIDS ARE RUNNING 1-WAY COUPLED.  IF SO, RETURN
      LGRCNVG = 1
      DO LG=2,NGRIDS
        CALL SGWF2LGR2PNT(LG)            
        IF(MXLGRITER .GT. 1)THEN
          LGRCNVG = 0
          EXIT 
        ENDIF
      ENDDO
      IF(LGRCNVG .EQ. 1) RETURN
C
C2------NOT RUNNING 1-WAY COUPLED SO LOOP THROUGH EACH GRID AND
C2------A: IF REQUESTED, WRITE CONVERGENCE OF LGR ITERATIONS  
C2------B: CHECK CONVERGENCE OF EACH GRID
C2------C: EXIT AND PRINT A MESSAGE IF MXLGRITER IS EXCEEDED
      LGRCNVG = 1
      IOFLG=0
      DO LG=2,NGRIDS
        CALL SGWF2LGR2PNT(LG)            
        CALL SGWF2BAS7PNT(LG)            
C2A
        IF(IOUTLGR .LT. 0 .AND. IOFLG .EQ. 0)THEN
          WRITE(*,500) LGRITER
          IOFLG=1
        ENDIF
        IF(IOUTLGR .LT. 0)THEN
          WRITE(*,501) LG
          WRITE(*,505) HCLOSELGR,HDIFFM,NODEH(1),NODEH(2),NODEH(3),
     &                 GNHEAD(NCMAXH,IBMAXH)
          WRITE(*,510) FCLOSELGR,FDIFFM,NODEF(1),NODEF(2),NODEF(3),
     &                 GNFLUX(NCMAXF,IBMAXF)
        ELSEIF(IOUTLGR .GT. 0)THEN
          WRITE(IOUT,500) LGRITER
          WRITE(IOUT,505) HCLOSELGR,HDIFFM,NODEH(1),NODEH(2),NODEH(3),
     &                    GNHEAD(NCMAXH,IBMAXH)
          WRITE(IOUT,510) FCLOSELGR,FDIFFM,NODEF(1),NODEF(2),NODEF(3),
     &                    GNFLUX(NCMAXF,IBMAXF)
        ENDIF
C2
C2B-----CHECK FOR CONVERGENCE IF NOT ON THE VERY FIRST LGR ITERATION 
C2B-----OF THE FIRST STRESS PERIOD OF THE FIRST TIME STEP. 
C2B-----(FORCE INITIAL LGR ITERATION)  
        IF(LGRITER .GT. 1 .OR. KPER .NE. 1 .OR. KSTP .NE. 1)THEN
          IF(ABS(HDIFFM) .GT. HCLOSELGR .OR. 
     &       ABS(FDIFFM) .GT. FCLOSELGR) LGRCNVG = 0
        ELSE
          LGRCNVG = 0
        ENDIF
!swm                          lgrcnvg=0   !swm:hack to force mxlgriter
C2 
C2C-----IF MXLGRITER IS EXCEEDED, END ITERATION AND PRINT MESSAGE
        IF(LGRITER .GE. MXLGRITER)THEN
          LGRCNVG = 2
          IF(IOUTLGR .LT. 0)WRITE(*,515) LG
          IF(IOUTLGR .GE. 0)WRITE(IOUT,515) LG
        ENDIF
      ENDDO
C3
      RETURN
      END SUBROUTINE GWF2LGR2CNVG
      
C***********************************************************************

*C-----VERSION 1.2 24JUNE2009 GWF2LGR2PBD
*      SUBROUTINE GWF2LGR2PBD(KSTP,KPER,IUBCF,IULPF,IUHUF,NGRIDS)
*C     ******************************************************************
*C     CALCULATE VOLUMETRIC BUDGET FOR PARENT B.C.  
*C     FOR COMPACT BUDGET FILES ALSO WRITE IFACE INFO
*C     NOTE: NEED TO HANDLE THE JOINT PARENT/CHILD CASE
*C     ******************************************************************
*C
*C        SPECIFICATIONS:
*C     ------------------------------------------------------------------
*      USE GLOBAL,      ONLY:IBOUND,IOUT,NCOL,NROW,NLAY
*      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,DELT,PERTIM,TOTIM,VBVL,VBNM
*      USE LGRMODULE,   ONLY:IUPBFSV,IUPBHSV,LGRDAT
*      USE GWFBCFMODULE,ONLY:IBCFCB
*      USE GWFLPFMODULE,ONLY:ILPFCB
*C      USE GWFHUFMODULE,ONLY:IHUFCB
*      CHARACTER(LEN=16):: TEXT, LGRAUX(1)
*      DOUBLE PRECISION RATIN,RATOUT,QQ,DZERO
*      DIMENSION GNFACE(1)
*      DATA TEXT /' GHOST-NODE FLUX'/
*      DATA LGRAUX /'           IFACE'/
*C     ------------------------------------------------------------------
*C
*C1------CLEAR RATIN AND RATOUT ACCUMULATORS. SET FLOW FILE UNIT NUMBER
*      DZERO=0.D0
*      RATIN=DZERO
*      RATOUT=DZERO
*      VBVL(3,MSUM)= 0.
*      VBVL(4,MSUM)= 0.
*      Q=0.
*      IF(IUBCF .NE. 0) IPLGRCB=IBCFCB
*      IF(IULPF .NE. 0) IPLGRCB=ILPFCB
*      IF(IUHUF .NE. 0) IPLGRCB=IHUFCB
*
*C2------LOOP THROUGH ALL SUBGRIDS
*      DO LG=2,NGRIDS
*C
*C3------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
*        IBD=ICBCFL
*        IF(IBD.EQ.2) THEN
*          NAUX=1   
*          CALL UBDSV4(KSTP,KPER,TEXT,NAUX,LGRAUX,IPLGRCB,NCOL,NROW,
*     1                NLAY,LGRDAT(LG)%NPBNODES,IOUT,DELT,PERTIM,TOTIM,
*     2                IBOUND)
*        END IF
*
*C4------REVERSE Q FOR PARENT PERSPECTIVE (+GNFLUX IS TO CHILD)
*C4------FLOW RATE IS POSITIVE (RECHARGE). ADD IT TO RATIN.
*C4------FLOW RATE IS NEGATIVE (DISCHARGE). ADD IT TO RATOUT.
*        QQ= SUM(-LGRDAT(LG)%GNFLUXOLD,-LGRDAT(LG)%GNFLUXOLD .GT. DZERO)
*        RATIN=RATIN+QQ
*        QQ= SUM(-LGRDAT(LG)%GNFLUXOLD,-LGRDAT(LG)%GNFLUXOLD .LT. DZERO)
*        RATOUT=RATOUT-QQ
*C            
*C4------CHECK IF WRITING COMPACT BUDGETS.  IF SO, LOOP THROUGH
*C4------PARENT CELLS AND SUM UP ALL GHOST-NODE FLUXES MATCHING THE
*C4------PARENT INDEX.  LOOP THROUGH ALL CHILD INTERFACE CELLS AND
*C4------FIND FIRST MATCH AND SWAP IFACE FOR PARENT PERSPECTIVE.
*C4------NOTE: THE COMPACT BUDGET NORMALLY SKIPS CELLS WITH 0 VALUES,
*C4------(ie INACTIVE CELLS) HOWEVER, WITH DRYING/WETTING THIS REQUIRES
*C4------KEEPING TRACK OF THE NUMBER OF ACTIVE CONNECTIONS.  FOR NOW, 
*C4------ALL INTERFACE CELLS ARE WRITTEN.  THIS WORKS BECAUSE GNFLUXOLD
*C4------IS ZEROED IN GWF2LGR2FMBF AND INACTIVE CELLS ARE SKIPPED. 
*        IF(IBD.EQ.2) THEN 
*          DO KP = 1, NLAY
*            DO IP = 1, NROW
*              DO JP = 1, NCOL
*                Q =SUM(-LGRDAT(LG)%GNFLUXOLD,LGRDAT(LG)%KPLC .EQ. KP 
*     1                       .AND. LGRDAT(LG)%IPLC .EQ. IP  
*     2                       .AND. LGRDAT(LG)%JPLC .EQ. JP)  
*                DO IB = 1,LGRDAT(LG)%NBNODES 
*                  IF(LGRDAT(LG)%KPLC(1,IB) .EQ. KP .AND. 
*     1               LGRDAT(LG)%IPLC(1,IB) .EQ. IP .AND.
*     2               LGRDAT(LG)%JPLC(1,IB) .EQ. JP) THEN
*
*                    IFACE = LGRDAT(LG)%IFACEGN(1,IB)
*                    IF (IFACE .EQ. 1) IFACE =2
*                    IF (IFACE .EQ. 2) IFACE =1
*                    IF (IFACE .EQ. 3) IFACE =4
*                    IF (IFACE .EQ. 4) IFACE =3
*                    IF (IFACE .EQ. 5) IFACE =6
*                    IF (IFACE .EQ. 6) IFACE =5
*                    GNFACE(1)=IFACE
*                    CALL UBDSVB(IPLGRCB,NCOL,NROW,JP,IP,KP,Q,
*     1                          GNFACE,1,NAUX,1,IBOUND,NLAY)
*                    Q=0.
*C-------MATCH FOUND. SKIP TO NEXT JP INDEX TO AVOID REPEATS
*                    EXIT 
*                  ENDIF
*                ENDDO
*              ENDDO
*            ENDDO
*          ENDDO
*        ENDIF 
*
*      ENDDO
*C
*C5------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
*      RIN=RATIN
*      ROUT=RATOUT
*      VBVL(3,MSUM)=VBVL(3,MSUM)+RIN
*      VBVL(4,MSUM)=VBVL(4,MSUM)+ROUT
*      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
*      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
*      VBNM(MSUM)=TEXT
*C
*C6------INCREMENT BUDGET TERM COUNTER(MSUM) 
*      MSUM=MSUM+1
*C
*C7------CALL ROUTINE TO OUTPUT PARENT B.C.'s FOR BFH IF REQUESTED
*      IF(IUPBHSV .NE. 0 .OR. IUPBFSV .NE. 0)
*     & CALL SGWF2LGR2BFHPOT(KSTP,KPER,NGRIDS)
*C8------RETURN
*      RETURN
*      END SUBROUTINE GWF2LGR2PBD
      
C***********************************************************************

C-----VERSION 1.2 22JUNE2009 SGWF2LGR2BFHPOT
      SUBROUTINE SGWF2LGR2BFHPOT(KSTP,KPER,NGRIDS)
C     ******************************************************************
C     OUTPUT PARENT GRID INTERFACE BOUNDARY FLUXES AND HEADS FOR BFH
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,NSTP,HNEW,GLOBALDAT
      USE GWFBASMODULE,ONLY:TOTIM
      USE LGRMODULE,   ONLY:ISCHILD,IUPBHSV,IUPBFSV,IBPFLG,LGRDAT
      CHARACTER(LEN=17):: TEXTH,TEXTF,TEXTR
      DATA TEXTH /'PARENT HEAD CHECK'/
      DATA TEXTF /'  GHOST-NODE FLUX'/
C     ------------------------------------------------------------------
C
C1------IF ON FIRST ITERATION, WRITE HEADER INFO AND NODE LOCATIONS FOR 
C1------INTERFACE CELLS
      IF(KPER .EQ. 1 .AND. KSTP .EQ. 1)THEN
        NPBNODES2 = 0
        NCBNODES2 = 0
        DO LG = 2, NGRIDS
          NPBNODES2 = NPBNODES2 + LGRDAT(LG)%NPBNODES
          NCBNODES2 = NCBNODES2 + SUM(LGRDAT(LG)%NCON)
        ENDDO
        TEXTR = 'SINGLE AREA' 
        IF(NGRIDS .GT. 2) TEXTR = 'MULTIPLE AREAS'
        IF(IUPBHSV .NE. 0)THEN 
          WRITE(IUPBHSV,200) TEXTR, NGRIDS
          WRITE(IUPBHSV,300) TEXTH,ISCHILD,NLAY,NROW,NCOL,SUM(NSTP),
     &                       NPBNODES2,NCBNODES2,IUPBHSV 
          WRITE(IUPBHSV,*) (IBPFLG(LG),LG=2,NGRIDS)
        ENDIF
        IF(IUPBFSV .NE. 0) THEN
          WRITE(IUPBFSV,200) TEXTR,NGRIDS
          WRITE(IUPBFSV,300) TEXTF,ISCHILD,NLAY,NROW,NCOL,SUM(NSTP),
     &                       NPBNODES2,NCBNODES2,IUPBHSV 
          WRITE(IUPBFSV,*) (IBPFLG(LG),LG=2,NGRIDS)
        ENDIF
C
        DO LG=2,NGRIDS
          IF (IUPBHSV .NE. 0) THEN
C2------LOOP THROUGH ALL PARENT CELLS AND WRITE HEADER FOR PARENT HEADS
            DO KP = 1,NLAY
              DO IP = 1,NROW
                DO JP = 1,NCOL
                  DO IB = 1, LGRDAT(LG)%NBNODES
                    IF(LGRDAT(LG)%KPLC(1,IB) .EQ. KP .AND. 
     1                 LGRDAT(LG)%IPLC(1,IB) .EQ. IP .AND.
     2                 LGRDAT(LG)%JPLC(1,IB) .EQ. JP) THEN
                       WRITE(IUPBHSV,400) KP,IP,JP,IBPFLG(LG)
C--------SKIP TO NEXT JP INDEX TO AVOID REPEATS
                      EXIT       
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO
            ENDDO
          ENDIF
C2A-----CHECK IF WRITING CHILD FLUXES.  IF SO, LOOP THROUGH ALL CHILD 
C2A-----CELLS AND WRITE HEADER FOR CHILD FLUXES.
          IF (IUPBFSV .NE. 0) THEN
            IB=0
            DO K=1,GLOBALDAT(LG)%NLAY
              DO I=1,GLOBALDAT(LG)%NROW
                DO J=1,GLOBALDAT(LG)%NCOL
                  IF (LGRDAT(LG)%ICBOUND(J,I,K) .EQ. IBPFLG(LG)) THEN
                    IB=IB+1            
                    DO NC=1,LGRDAT(LG)%NCON(IB)
C-------WRITE PARENT CELL INDICES FOLLOWED BY CHILD CELL INDICES  
C-------SWITCH IFACE TO BE FROM PARENT PERSPECTIVE                      
                      IFACE = LGRDAT(LG)%IFACEGN(NC,IB)
                      IF (IFACE .EQ. 1) IFACE =2
                      IF (IFACE .EQ. 2) IFACE =1
                      IF (IFACE .EQ. 3) IFACE =4
                      IF (IFACE .EQ. 4) IFACE =3
                      IF (IFACE .EQ. 5) IFACE =6
                      IF (IFACE .EQ. 6) IFACE =5
                      WRITE(IUPBFSV,450) LGRDAT(LG)%KPLC(NC,IB),
     1                                   LGRDAT(LG)%IPLC(NC,IB),
     2                                   LGRDAT(LG)%JPLC(NC,IB),
     3                                   K,I,J,IFACE,IBPFLG(LG)
                    ENDDO
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF   
C
C3------WRITE HEADER TIME STEP INFO
      IF(IUPBHSV .NE. 0) WRITE(IUPBHSV,500) KPER,KSTP,TOTIM
      IF(IUPBFSV .NE. 0) WRITE(IUPBFSV,500) KPER,KSTP,TOTIM
C
C4------LOOP THROUGH ALL GRIDS.  IF A BOUNDARY INTERFACE CELL THEN 
c4------WRITE BOUNDARY HEADS AND FLUXES IF DESIRED.  
      DO LG=2,NGRIDS
        IF (IUPBHSV .NE. 0) THEN
          DO KP = 1,NLAY
            DO IP = 1,NROW
              DO JP = 1,NCOL
                DO IB = 1, LGRDAT(LG)%NBNODES
                  IF(LGRDAT(LG)%KPLC(1,IB) .EQ. KP .AND. 
     1             LGRDAT(LG)%IPLC(1,IB) .EQ. IP .AND.
     2             LGRDAT(LG)%JPLC(1,IB) .EQ. JP) THEN
                    WRITE(IUPBHSV,600) HNEW(JP,IP,KP)
C--------SKIP TO NEXT JP INDEX TO AVOID REPEATS
                    EXIT       
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF 
C4A-----CHECK IF WRITING CHILD FLUXES.  IF SO LOOP THROUGH ALL CHILD 
C4A-----CELLS AND WRITE CHILD FLUXES
        IF (IUPBFSV .NE. 0) THEN
          IB=0
          DO K=1,GLOBALDAT(LG)%NLAY
            DO I=1,GLOBALDAT(LG)%NROW
              DO J=1,GLOBALDAT(LG)%NCOL
                IF (LGRDAT(LG)%ICBOUND(J,I,K) .EQ. IBPFLG(LG)) THEN
                  IB=IB+1            
                  DO NC=1,LGRDAT(LG)%NCON(IB)
C-------SWITCH SIGN:  GNFLUX IS DEFINED FROM THE CHILD PERSPECTIVE
C-------GNFLUXOLD WAS THE FLUX APPLIED TO THE PARENT THIS ITERATION
                    WRITE(IUPBFSV,600) -LGRDAT(LG)%GNFLUXOLD(NC,IB)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
C
 200  FORMAT(1X,A,1X,I4)
 300  FORMAT(1X,A,8(I4,2X))
 400  FORMAT(1X,4(I4,2X)) 
 450  FORMAT(1X,8(I4,2X)) 
 500  FORMAT(1X,'KPER=',I4,2x,'KSTP=',I4,2X,'TOTIM=',G14.7) 
 600  FORMAT(1X,G16.9) 
C
C5----RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2BFHPOT
      
C***********************************************************************

*C-----VERSION 1.2 24JUNE2009 GWF2LGR2CBD
*      SUBROUTINE GWF2LGR2CBD(KSTP,KPER,IUBCF,IULPF,IUHUF)
*C     ******************************************************************
*C     ACCUMULATE THE GHOST-NODE FLUXES ACROSS THE CHILD GRID
*C     FOR COMPACT BUDGET FILES ALSO WRITE IFACE INFO
*C     STORE PARENT AND CHILD GHOST-NODE FLUX RATES FOR PRINTING
*C     ******************************************************************
*C
*C        SPECIFICATIONS:
*C     ------------------------------------------------------------------
*      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IOUT,IBOUND,HNEW,GLOBALDAT
*      USE GWFBASMODULE,ONLY:MSUM,ICBCFL,DELT,PERTIM,TOTIM,VBVL,VBNM
*      USE LGRMODULE,   ONLY:PRATIN,CRATIN,PRATOUT,CRATOUT,IBFLG,NCON,
*     1                      IUCBHSV,IUCBFSV,MXLGRITER,JPLC,IPLC,KPLC,
*     2                      IFACEGN,GNFLUXOLD,GNFLUX,ICBOUND
*      USE GWFBCFMODULE,ONLY:IBCFCB
*      USE GWFLPFMODULE,ONLY:ILPFCB
**      USE GWFHUFMODULE,ONLY:IHUFCB
*      CHARACTER(LEN=16):: TEXT, LGRAUX(1)
*      DOUBLE PRECISION RATIN,RATOUT,PPRATIN,PPRATOUT,QQ,DZERO
*      DIMENSION GNFACE(1)
*      DATA TEXT /' GHOST-NODE HEAD'/
*      DATA LGRAUX /'           IFACE'/
*C     ------------------------------------------------------------------
*C
*C1------IF RUNNING 1-WAY COUPLED, DO NOT OUTPUT BOUNDARY FLUXES
*cswm: keep these so we get correct budget on child when 1-way coupled      IF(MXLGRITER .LE. 1) RETURN
*C      
*C2------CLEAR RATIN AND RATOUT ACCUMULATORS. SET FLOW FILE UNIT NUMBER
*      DZERO=0.D0
*      RATIN=DZERO
*      RATOUT=DZERO
*      PPRATIN=DZERO
*      PPRATOUT=DZERO
*      VBVL(3,MSUM)= 0.
*      VBVL(4,MSUM)= 0.
*      IF(IUBCF .NE. 0) ICLGRCB=IBCFCB
*      IF(IULPF .NE. 0) ICLGRCB=ILPFCB
*      IF(IUHUF .NE. 0) ICLGRCB=IHUFCB
*
*C3------IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
*        IBD=ICBCFL
*        IF(IBD.EQ.2) THEN
*          NCBFACES = SUM(NCON)
*          NAUX=1   
*          CALL UBDSV4(KSTP,KPER,TEXT,NAUX,LGRAUX,ICLGRCB,NCOL,NROW,
*     1                NLAY,NCBFACES,IOUT,DELT,PERTIM,TOTIM,IBOUND)
*        END IF
*
*C3------LOOP THROUGH ALL CHILD CELLS AND ACCUMULATE GHOST-NODE FLUXES
*C3------SUM UP ALL POSITIVE FLOW INTO RATIN (RECHARGE TO CHILD)
*C3------SUM UP ALL NEGATIVE FLOW INTO RATOUT (DISCHARGE FROM CHILD)
*C3------NOTE: THE COMPACT BUDGET NORMALLY SKIPS CELLS WITH 0 VALUES,
*C3------(ie INACTIVE CELLS) HOWEVER, WITH DRYING/WETTING THIS REQUIRES
*C3------KEEPING TRACK OF THE NUMBER OF ACTIVE CONNECTIONS.  FOR NOW, 
*C3------ALL INTERFACE CELLS ARE WRITTEN.  THIS WORKS BECAUSE GNFLUX
*C3------IS ZEROED IN GWF2LGR2FMBF AND INACTIVE CELLS ARE SKIPPED. 
*      IB = 0
*      DO K = 1, NLAY
*        DO I = 1, NROW
*          DO J = 1, NCOL
*            IF(ICBOUND(J,I,K) .EQ. IBFLG) THEN
*              IB=IB+1
*C-------LOOP THROUGH ALL CONNECTIONS AND ACCUMULATE FLUXES
*              DO NC=1,NCON(IB)
*                Q=GNFLUX(NC,IB)
*                QQ=Q
*                IF(Q .GT. 0.) RATIN=RATIN+QQ
*                IF(Q .LT. 0.) RATOUT=RATOUT-QQ
*C
*C3A-----CHECK IF WRITING COMPACT BUDGETS.  IF SO, SET IFACE
*C3A-----AND WRITE GHOST-NODE FLUXES AND GNFACE
*                IF(IBD.EQ.2) THEN 
*                  GNFACE(1)=IFACEGN(NC,IB)
*                  CALL UBDSVB(ICLGRCB,NCOL,NROW,J,I,K,Q,GNFACE,1,NAUX,
*     1                        1,IBOUND,NLAY)
*                ENDIF
*C
*C3B-----SUM UP FLUXES FROM PARENT FOR COMPARISON
*                Q=GNFLUXOLD(NC,IB)     
*                QQ=Q
*                IF(Q .GT. 0.) PPRATIN=PPRATIN+QQ
*                IF(Q .LT. 0.) PPRATOUT=PPRATOUT-QQ
*C
*              ENDDO
*            ENDIF
*          ENDDO
*        ENDDO
*      ENDDO
*C
*C4------MOVE RATES, VOLUMES & LABELS INTO ARRAYS FOR PRINTING.
*      RIN=RATIN
*      ROUT=RATOUT
*      CRATIN=RATIN
*      CRATOUT=RATOUT
*      PRATIN=PPRATIN
*      PRATOUT=PPRATOUT
*      VBVL(3,MSUM)=VBVL(3,MSUM)+RIN
*      VBVL(4,MSUM)=VBVL(4,MSUM)+ROUT
*      VBVL(1,MSUM)=VBVL(1,MSUM)+RIN*DELT
*      VBVL(2,MSUM)=VBVL(2,MSUM)+ROUT*DELT
*      VBNM(MSUM)=TEXT
*C
*C5------INCREMENT BUDGET TERM COUNTER(MSUM) 
*      MSUM=MSUM+1
*C
*C6------CALL ROUTINE TO OUTPUT CHILD B.C.'s FOR BFH IF REQUESTED
*      IF(IUCBHSV .NE. 0 .OR. IUCBFSV .NE. 0) CALL SGWF2LGR1BFHCOT(KSTP,
*     &  KPER,GLOBALDAT(1)%NCOL,GLOBALDAT(1)%NROW,GLOBALDAT(1)%NLAY)
*C7----RETURN
*      RETURN
*      END SUBROUTINE GWF2LGR2CBD
      
C***********************************************************************

C-----VERSION 1.0 16JULY2010 GWF2LGR2COT
      SUBROUTINE GWF2LGR2COT(KSTP,KPER,IGRID)
C     ******************************************************************
C     PRINT PARENT AND CHILD GRID GHOST-NODE FLUXES FOR COMPARISON
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IOUT
      USE GWFBASMODULE,ONLY:TOTIM
      USE LGRMODULE,   ONLY:MXLGRITER,IOUTLGR,PRATIN,CRATIN,PRATOUT,
     &                      CRATOUT
C     ------------------------------------------------------------------
C      
C1------IF RUNNING 1-WAY COUPLED, DO NOT OUTPUT GHOST-NODE FLUXES
      IF(MXLGRITER .LE. 1) RETURN
C
C2------INITIALIZE VARIABLES
      PDIFFRIN=0.
      PDIFFROUT=0.
C3------CALCULATE DIFFERENCE BETWEEN PARENT AND CHILD GHOST-NODE FLUXES
      DIFFRIN = PRATIN-CRATIN
      AVGRAT = (PRATIN+CRATIN)/2.
      IF(AVGRAT .NE. 0.) PDIFFRIN=100.*DIFFRIN/AVGRAT
      DIFFROUT = PRATOUT-CRATOUT
      AVGRAT = (PRATOUT+CRATOUT)/2.
      IF(AVGRAT .NE. 0.) PDIFFROUT=100.*DIFFROUT/AVGRAT
C
C4------PRINT DIFFERENCES TO OUTPUT FILE.  
C
      WRITE(IOUT,200) KSTP, KPER
      WRITE(IOUT,210)
      WRITE(IOUT,220) PRATIN, CRATIN, DIFFRIN, PDIFFRIN
      WRITE(IOUT,230) PRATOUT, CRATOUT, DIFFROUT, PDIFFROUT
C5------PRINT TO SCREEN IF REQUESTED 
      IF(IOUTLGR .LT. 0) THEN
        WRITE(*,205) IGRID
        WRITE(*,210)
        WRITE(*,220) PRATIN, CRATIN, DIFFRIN, PDIFFRIN
        WRITE(*,230) PRATOUT, CRATOUT, DIFFROUT, PDIFFROUT
      END IF
C
  200 FORMAT('1',/2X,'FLUX ACROSS PARENT-CHILD INTERFACE AT'
     1,' TIME STEP',I3,' IN STRESS PERIOD ',I4/2X,72('-'))
  205 FORMAT(/,1X,'GRID NUMBER=',I3)
  210 FORMAT(1X,'G-N FLUX ',7X,'PARENT ',6X,' CHILD ', 8X,'DIFFERENCE',
     &       6X,'% DIFFERENCE')
  220 FORMAT(1X,' RATE IN:'4X,2(ES12.4,2X),1X, ES12.4, 2X,F12.4)
  230 FORMAT(1X,'RATE OUT:'4X,2(ES12.4,2X),1X, ES12.4, 2X,F12.4)
C
C6----RETURN
      RETURN
      END SUBROUTINE GWF2LGR2COT
C
C***********************************************************************

C-----VERSION 1.0 22JUNE2009 SGWF2LGR1BFHCOT
      SUBROUTINE SGWF2LGR1BFHCOT(KSTP,KPER,NCOLP,NROWP,NLAYP)
C     ******************************************************************
C     OUTPUT CHILD GRID INTERFACE BOUNDARY FLUXES AND HEADS FOR BFH 
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,NSTP
      USE GWFBASMODULE,ONLY:TOTIM
      USE LGRMODULE,   ONLY:ISCHILD,NCPP,IBOTFLG,IBFLG,IUCBHSV,IUCBFSV,
     &                      GNFLUX,GNHEAD,GNCOND,IFACEGN,KPLC,IPLC,JPLC,
     &                      NCON,ICBOUND,NBNODES,NPBNODES
      CHARACTER(LEN=17):: TEXTH,TEXTF
      DATA TEXTH /'  GHOST-NODE HEAD'/
      DATA TEXTF /' CHILD FLUX CHECK'/
C     ------------------------------------------------------------------
C
C1------IF ON FIRST ITERATION, WRITE HEADER INFO, AND NODE LOCATIONS 
C1------FOR INTERFACE CELLS.  NOTE THAT HEADS ARE ASSOCIATED WITH CHILD 
C1------GRID INDICES WHILE FLUXES ARE ASSOCIATED WITH PARENT GRID INDICES
      IF(KPER .EQ. 1 .AND. KSTP .EQ. 1)THEN
C2------FIND THE TOTAL NUMBER OF GHOST NODE CONNECTIONS 
          NBFACES = SUM(NCON)
C2------CHECK IF SAVING BOUNDARY HEADS.
        IF(IUCBHSV .NE. 0)THEN 
          WRITE(IUCBHSV,300) TEXTH,ISCHILD,NLAY,NROW,NCOL,SUM(NSTP),
     &                      NPBNODES,NBFACES,NCPP,IBOTFLG,IBFLG,IUCBFSV
          IB=0
          DO K=1,NLAY
            DO I=1,NROW
              DO J=1,NCOL
                IF(ICBOUND(J,I,K) .EQ. IBFLG)THEN 
                  IB=IB+1
                  DO NC=1,NCON(IB)
                    WRITE(IUCBHSV,400) K,I,J,IFACEGN(NC,IB),KPLC(NC,IB),
     &                                 IPLC(NC,IB),JPLC(NC,IB),IBFLG  
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF   
C
C3------CHECK IF SAVING BOUNDARY FLOWS.
C3------WRITE CELL INDICES FOR CHILD GRID FLUXES AND THEIR CORESPONDING 
C3------LOCATIONS IN THE PARENT GRID.
        IF(IUCBFSV .NE. 0)THEN
          WRITE(IUCBFSV,300) TEXTF,ISCHILD,NLAYP,NROWP,NCOLP,SUM(NSTP),
     &                      NPBNODES,NBFACES,NCPP,IBOTFLG,IBFLG,IUCBFSV
          IB=0
          DO K=1,NLAY
            DO I=1,NROW
              DO J=1,NCOL
                IF(ICBOUND(J,I,K) .EQ. IBFLG) THEN
                  IB=IB+1
                  DO NC = 1,NCON(IB)
                    WRITE(IUCBFSV,400) K,I,J,IFACEGN(NC,IB),KPLC(NC,IB),
     &                                 IPLC(NC,IB),JPLC(NC,IB),IBFLG  
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF   
C 
      ENDIF   !END FIRST TIME STEP OF FIRST STRESS PERIOD 
C
C4------IF BOUNDARY HEADS OR FLUXES ARE BEING SAVED, WRITE 
C4------CORRESPONDING HEADER TIME STEP INFO.
C4------LOOP THROUGH ALL CELLS.  IF A BOUNDARY INTERFACE CELL THEN 
C4------WRITE BOUNDARY HEAD AND/OR FLUXES, IF REQUESTED 
      IF(IUCBHSV .NE. 0 .OR. IUCBFSV .NE. 0)THEN
        IF(IUCBHSV .NE. 0) WRITE(IUCBHSV,500) KPER,KSTP,TOTIM 
        IF(IUCBFSV .NE. 0) WRITE(IUCBFSV,500) KPER,KSTP,TOTIM 
        IB=0
        DO K=1,NLAY
          DO I=1,NROW
            DO J=1,NCOL
              IF(ICBOUND(J,I,K) .EQ. IBFLG)THEN 
                IB=IB+1
                DO NC=1,NCON(IB)
                  IF(IUCBHSV .NE. 0) WRITE(IUCBHSV,600) GNHEAD(NC,IB), 
     1                                                  GNCOND(NC,IB)
                  IF(IUCBFSV .NE. 0) WRITE(IUCBFSV,600) GNFLUX(NC,IB)
                ENDDO
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
C
 300  FORMAT(1X,A,5(I4,2X),2(I6,2x),4(I4,2X)) 
 400  FORMAT(1X,8(I4,2X)) 
 500  FORMAT(1X,'KPER=',I4,2x,'KSTP=',I4,2X,'TOTIM=',G14.7) 
 600  FORMAT(1X,2G16.9) 
C
C6----RETURN
      RETURN
      END SUBROUTINE SGWF2LGR1BFHCOT
      
C***********************************************************************

C-----VERSION 1.2 06JANUARY2009 SGWF2LGR2GNCOND
      SUBROUTINE SGWF2LGR2GNCOND(I,J,K,JP,IP,KP,IUBCF,IULPF,IUHUF,
     1                           IUBCFC,IULPFC,IUHUFC,IB,NC,LAYHDTFLG,
     2                           THICKP,THICKPC,DELRP,DELCP,THICKC,LG,
     3                           CGN)
C     ******************************************************************
C     CALCULATE THE GHOST-NODE CONDUCTANCE OF THE MATERIAL BETWEEN A
C     PARENT AND A CHILD CELL.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:GLOBALDAT
      USE LGRMODULE,   ONLY:LGRDAT
C     ------------------------------------------------------------------

C1------SET PARENT HK's  (NOTE: hardwired for parent grid)
      CALL SGWF2LGR2HKP(JP,IP,KP,IUBCF,IULPF,IUHUF,HKPX,HKPY,HKPZ,1)
C2------SET CHILD HK'S 
      CALL SGWF2LGR2HKC(J,I,K,IUBCFC,IULPFC,IUHUFC,HKX,HKY,HKZ,LG) 
C
C3------IF PARENT OR CHILD CELL IS DEWATERED, THEN SET CGN = 0 AND
C3------RETURN
      IF (THICKP .LE. 0. .OR. THICKC .LE. 0.) THEN
        CGN = 0.
        RETURN
      END IF
C
C4------CHECK CONNECTION FACE TO DETERMINE ORIENTATION
      IF (LGRDAT(LG)%IFACEGN(NC,IB) .EQ. 1 .OR.
     1  LGRDAT(LG)%IFACEGN(NC,IB) .EQ. 2 )THEN
        WIDTH=GLOBALDAT(LG)%DELC(I)
C-------CHECK IF PARENT IS UNCONFINED TO USE APPROPRIATE HEIGHT
        IF (LAYHDTFLG .EQ. 0) THEN
          HEIGHTP=THICKP/LGRDAT(LG)%NCPPL(KP) 
        ELSE
          IF (THICKPC .GT. 0) THEN
            HEIGHTP=THICKPC
          ELSE
            HEIGHTP=0.
          ENDIF
        END IF
        HEIGHTC=THICKC
        DP = 0.5*DELRP
        DC = 0.5*GLOBALDAT(LG)%DELR(J)
        HKP = HKPX
        HKC = HKX

      ELSE IF (LGRDAT(LG)%IFACEGN(NC,IB) .EQ. 3 .OR.
     1         LGRDAT(LG)%IFACEGN(NC,IB) .EQ. 4 )THEN
        WIDTH=GLOBALDAT(LG)%DELR(J)
        IF (LAYHDTFLG .EQ. 0) THEN
          HEIGHTP=THICKP/LGRDAT(LG)%NCPPL(KP) 
        ELSE
          IF (THICKPC .GT. 0) THEN
            HEIGHTP=THICKPC
          ELSE
            HEIGHTP=0.
          ENDIF
        END IF
        HEIGHTC=THICKC
        DP = 0.5*DELCP
        DC = 0.5*GLOBALDAT(LG)%DELC(I)
        HKP = HKPY
        HKC = HKY
      
      ELSE IF (LGRDAT(LG)%IFACEGN(NC,IB) .EQ. 5 .OR.
     1         LGRDAT(LG)%IFACEGN(NC,IB) .EQ. 6 )THEN
        WIDTH=GLOBALDAT(LG)%DELR(J)
        HEIGHTP=GLOBALDAT(LG)%DELC(I)
        HEIGHTC=GLOBALDAT(LG)%DELC(I)
        DP = 0.5*THICKP
        DC = 0.5*THICKC
        HKP = HKPZ
        HKC = HKZ
      END IF
C5------FORMULATE CONDUCTANCE FROM PARENT SIDE AND CHILD SIDE.
C5------USE THESE CONDUCTANCES IN SERIES TO GET GHOST-NODE CONDUCTANCE
      CONDP = HKP*WIDTH*HEIGHTP/DP
      CONDC = HKC*WIDTH*HEIGHTC/DC
      CGN = CONDP*CONDC/(CONDP+CONDC)
       
!swmdbg       if (i .eq. 1 .and. j .eq. 1)                                   
!swmdbg     &         write(73,*) ib, k,i,j,nc,kp,ip,jp,                    
!swmdbg     1                     LGRDAT(LG)%IFACEGN(NC,IB),width,heightp,  
!swmdbg     2         heightc, dp,dc, hkp,hkc,condp, condc, cgn             
C
C6------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2GNCOND
      
C***********************************************************************

C-----VERSION 1.2 20NOVEMBER2007 SGWF2LGR12COL
      SUBROUTINE SGWF2LGR2PLOC(NCPP,LBEG,LEND,L,NPLOC)
C     ******************************************************************
C     FIND CELL LOCATION WHERE PARENT SHARES AN INTERFACE WITH THE CHILD
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
C1------INITIALIZE CHILD CELL LOCATIONS
      LCBEG=1
      LCEND=LCBEG+NCPP-1
C2------LOOP THROUGH PARENT CELLS.  CHECK IF CURRENT CELL INDEX (L)
C2------FALLS WITHIN THIS PARENT CELL.  IF SO, STORE THE LOCATION AND 
C2------EXIT.  IF NOT, SHIFT TO NEXT PARENT CELL AND CONTINUE LOOP.
      
      DO LP=LBEG,LEND
        IF(L .GE. LCBEG .AND. L .LE. LCEND)THEN
          NPLOC=LP
          EXIT    
        ENDIF
        LCBEG=LCEND+1
        LCEND=LCBEG+NCPP-1
      ENDDO
C3------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2PLOC
      
C***********************************************************************

C-----VERSION 1.2 16JANUARY008 SGWF2LGR2HKP
      SUBROUTINE SGWF2LGR2HKP(JP,IP,KP,IUBCF,IULPF,IUHUF,HKPX,HKPY,
     1                        HKPZ,LG)
C     ******************************************************************
C     FIND PARENT X, Y, AND Z DIRECTION HYDRAULIC CONDUCTIVITIES 
C     DEPENDING ON WHICH FLOW PACKAGE IS ACTIVE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------

C1------CHECK WHICH FLOW PACKAGE IS USED AND CALL CORRESPONDING ROUTINE
      IF (IUBCF .NE. 0) CALL SGWF2LGR2BCFHK(JP,IP,KP,HKPX,HKPY,HKPZ,LG)
      IF (IULPF .NE. 0) CALL SGWF2LGR2LPFHK(JP,IP,KP,HKPX,HKPY,HKPZ,LG)
      !IF (IUHUF .NE. 0) CALL SGWF2LGR2HUFHK(JP,IP,KP,HKPX,HKPY,HKPZ,LG)
C2------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2HKP
      
C***********************************************************************

C-----VERSION 1.2 16JANUARY008 SGWF2LGR2HKC
      SUBROUTINE SGWF2LGR2HKC(J,I,K,IUBCF,IULPF,IUHUF,HKX,HKY,HKZ,LG)
C     ******************************************************************
C     FIND CHILD X, Y, AND Z DIRECTION HYDRAULIC CONDUCTIVITIES 
C     DEPENDING ON WHICH FLOW PACKAGE IS ACTIVE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------

C1------CHECK WHICH FLOW PACKAGE IS USED AND CALL CORRESPONDING ROUTINE
      IF (IUBCF .NE. 0) CALL SGWF2LGR2BCFHK(J,I,K,HKX,HKY,HKZ,LG)
      IF (IULPF .NE. 0) CALL SGWF2LGR2LPFHK(J,I,K,HKX,HKY,HKZ,LG)
      !IF (IUHUF .NE. 0) CALL SGWF2LGR2HUFHK(J,I,K,HKX,HKY,HKZ,LG)
C2------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2HKC
      
C***********************************************************************
C***********************************************************************

C-----VERSION 1.2 21JULY2010 SGWF2LGR2BCFHK
      SUBROUTINE SGWF2LGR2BCFHK(J,I,K,HKX,HKY,HKZ,IGRID)
C     ******************************************************************
C     FIND X, Y, AND Z DIRECTION HYDRAULIC CONDUCTIVITIES FROM THE BCF 
C     PACKAGE  
C     NOTE: VK NOT DEFINED.  APPROXIMATED BASED LEAKANCE AND VERTICAL
C     DISTANCE.  VALUE IS AVERAGE VK BETWEEN CELLS, NOT VK OF THE CELL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFBCFMODULE, ONLY:TRPY  !swm: NEED TO ADD TO BCF
      USE LGRMODULE,    ONLY:LGRDAT
      use GWFBCFMODULE, only: SGWF2BCF7PNT
C     ------------------------------------------------------------------
      CALL SGWF2BCF7PNT(IGRID)
C1------CALCULATE HORIZONTAL HYDRAULIC CONDUCTIVITIES.  HKY DEPENDS ON 
C1------ANISOTROPY
      HKX=LGRDAT(IGRID)%HK(J,I,K)
      HKY=HKX*TRPY(K)
C2------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY
      HKZ = LGRDAT(IGRID)%VK(J,I,K)      
C3------RETURN      
      RETURN
      END SUBROUTINE SGWF2LGR2BCFHK
      
C***********************************************************************

C-----VERSION 1.2 16JANUARY008 SGWF2LGR2LPFHK
      SUBROUTINE SGWF2LGR2LPFHK(J,I,K,HKX,HKY,HKZ,IGRID)
C     ******************************************************************
C     FIND X, Y, AND Z DIRECTION HYDRAULIC CONDUCTIVITIES FROM THE LPF 
C     PACKAGE  
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFLPFMODULE, ONLY:HK,CHANI,HANI,LAYVKA,VKA, 
     &                       SGWF2LPF7PNT
C     ------------------------------------------------------------------
      CALL SGWF2LPF7PNT(IGRID)
C1------CALCULATE HORIZONTAL HYDRAULIC CONDUCTIVITIES.  HKY DEPENDS ON 
C1------ANISOTROPY
      HKX=HK(J,I,K)
      IF (CHANI(K) .GT. 0.) THEN
        HKY=HKX*CHANI(K)
      ELSE
        HKY=HKX*HANI(J,I,K)
      END IF
C
C2------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY
      IF(LAYVKA(K).EQ.0) THEN
        HKZ=VKA(J,I,K)
      ELSE
        HKZ=HK(J,I,K)/VKA(J,I,K)
      END IF
C
C3------RETURN      
      RETURN
      END SUBROUTINE SGWF2LGR2LPFHK
      
C***********************************************************************

*C-----VERSION 1.2 22JULY2010 SGWF2LGR2HUFHK
*      SUBROUTINE SGWF2LGR2HUFHK(J,I,K,HKX,HKY,HKZ,IGRID)
*C     ******************************************************************
*C     FIND X, Y, AND Z DIRECTION HYDRAULIC CONDUCTIVITIES FROM THE HUF 
*C     PACKAGE  
*C     ******************************************************************
*C
*C        SPECIFICATIONS:
*C     ------------------------------------------------------------------
*      USE GLOBAL,       ONLY:NCOL,NROW,NLAY,IOUT,NBOTM,HNEW,LBOTM,BOTM,
*     1                       IBOUND
**      USE GWFHUFMODULE, ONLY:NHUF,LTHUF,HGUVANI,HUFHK,HUFVK,HUFKDEP,GS,
**     1                       VKAH,HK,HKCC,HUFTHK 
*      USE LGRMODULE,    ONLY:LGRDAT
*      REAL    MULTKDEP
*C     ------------------------------------------------------------------
*      CALL SGWF2HUF7PNT(IGRID)
*C1------CALCULATE HORIZONTAL HYDRAULIC CONDUCTIVITIES.  
*      HKX=HK(J,I,K)
*      HKY=HKCC(J,I,K)
*C
*C2------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY
*C2A-----CHECK IF CONFINED.  IF SO, USE VALUE STORED IN VK
*      IF(LTHUF(K) .EQ. 0) THEN
*        HKZ=LGRDAT(IGRID)%VK(J,I,K)        
*      ELSE
*C2B-----UNCONFINED.  CALCULATE UPDATED VALUE BASED ON SGWFHUF7VKL
*        DO 100 NU=1,NHUF
*C         FIND TOP AND BOTTOM LAYERS THIS UNIT APPLIES TO
*          TOPU=HUFTHK(J,I,NU,1)
*          THCKU=HUFTHK(J,I,NU,2)
*cswm          IF(ABS(THCKU).LT.1E-4) GOTO 210
*          IF(THCKU.EQ. 0.) CYCLE
*          BOTU=TOPU-THCKU
*C-----------Determine which layer(s) unit applies to
*          IFLG=1
*          CALL SGWF2HUF7HSRCH(NCOL,NROW,NLAY,BOTM,NBOTM,I,J,TOPU,
*     &                        BOTU,HNEW,IBOUND,KT,KB,IFLG)
*          IF(IFLG.EQ.1) CYCLE
*          IF(KT.GT.K) CYCLE
*          IF(KB.LT.K) CYCLE
*C         IF IN THIS LAYER, CONTRIBUTE TO VK
*C         GET THICKNESS OF UNIT IN THIS LAYER
*          CALL UHUF7THK(BOTM(J,I,LBOTM(K)-1),BOTM(J,I,LBOTM(K)),
*     &                      TOPU,THCKU,THCK,ATPU,ABTU)
*          IF(THCK.EQ.0.0) CYCLE
*          HUFHK(NU)=0.0
*          HUFVK(NU)=0.0
*          HUFKDEP(NU)=0.0
*          IF(HGUVANI(NU).EQ.0.) THEN
*            CALL UHUF7POP(HUFVK,'VK  ',I,J,NU,IOUT)
*            LGRDAT(IGRID)%VK(J,I,K)=LGRDAT(IGRID)%VK(J,I,K)+
*     &                              THCK/HUFVK(NU)
*          ELSE
*            MULTKDEP = 1.0
*            CALL UHUF7POP(HUFHK,'HK  ',I,J,NU,IOUT)
*            CALL UHUF7POP(HUFVK,'VANI',I,J,NU,IOUT)
*            CALL UHUF7POP(HUFKDEP,'KDEP',I,J,NU,IOUT)
*            IF(HUFKDEP(NU).NE.0.) CALL SGWF2HUF7KDEP(HUFKDEP(NU),TOPU,
*     &                                               BOTU,GS(J,I),
*     &                                               MULTKDEP)
*            LGRDAT(IGRID)%VK(J,I,K)=LGRDAT(IGRID)%VK(J,I,K)+
*     &                              THCK*HUFVK(NU)/(MULTKDEP*HUFHK(NU))
*          END IF
*  100   CONTINUE
*        HKZ=LGRDAT(IGRID)%VK(J,I,K)       
*      END IF
*C
*C3------RETURN      
*      RETURN
*      END SUBROUTINE SGWF2LGR2HUFHK
      
C***********************************************************************

      SUBROUTINE SGWF2LGR2XYZLOC(J,I,K,NCOL,NROW,NLAY,DELR,DELC,BOTM,X,
     1                           Y,Z)
C     ******************************************************************
C     FIND X, Y, AND Z COORDINATES OF CELL CENTERS.
C     (Y MEASURED FROM UPPER LEFT AND Z MEASURED FROM BOTTOM)
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION DELR(NCOL), DELC(NROW), BOTM(NCOL,NROW,0:NLAY)
C     ------------------------------------------------------------------
      X = SUM(DELR(1:J)) - 0.5*DELR(J)
      Y = SUM(DELC(1:I)) - 0.5*DELC(I)
      Z = BOTM(J,I,K) + 0.5*(BOTM(J,I,K-1)-BOTM(J,I,K))
      RETURN
      END SUBROUTINE SGWF2LGR2XYZLOC
      
C***********************************************************************

C-----VERSION 2.0 27JULY2012 SGWF2LGR2INTERIOR
      SUBROUTINE SGWF2LGR2INTERIOR(IGRID,NGRIDS,J,I,LINT)
C     ******************************************************************
C     DETERMINE IF THE CELL INDEX IS INSIDE OF THE PARENT AREA THAT IS
C     COVERED BY A CHILD.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE LGRMODULE,   ONLY:LGRDAT
      INTEGER LINT
C     ------------------------------------------------------------------
C1------INITIALIZE LINT TO ZERO.  
      LINT = 0
C
C2------SEARCH THROUGH ALL SUBGRIDS AND DETERMINE IF CELL IS INSIDE AN 
C2------AREA COVERED BY A CHILD.  IF SO, SET LINT=1 AND EXIT SEARCH. 
C NOTE: WILL NOT WORK FOR IRREGULARLY SHAPED AREAS
      DO LG = IGRID+1,NGRIDS
        IF (I .GE. LGRDAT(LG)%NPRBEG .AND. I .LE. LGRDAT(LG)%NPREND 
     &     .AND. 
     &     J .GE. LGRDAT(LG)%NPCBEG .AND. J .LE. LGRDAT(LG)%NPCEND) THEN
          LINT = 1
          EXIT
        ENDIF
      END DO
C
C3------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2INTERIOR
      
C***********************************************************************

C-----VERSION 2.0 15JUNE2013 SGWF2LGR2EVT
      SUBROUTINE SGWF2LGR2EVT(J,I)
C     ******************************************************************
C     ZERO OUT EVAPOTRANSPIRATION IN THE PARENT AREA THAT IS COVERED BY 
C     A CHILD.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFEVTMODULE,   ONLY:EVTR
C     ------------------------------------------------------------------
C1------ZERO OUT EVAPOTRANSPIRATION 
      EVTR(J,I) = 0.
C2------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2EVT
      
C***********************************************************************

C-----VERSION 2.0 15JUNE2013 SGWF2LGR2RCH
      SUBROUTINE SGWF2LGR2RCH(J,I)
C     ******************************************************************
C     ZERO OUT RECHARGE IN THE PARENT AREA THAT IS COVERED BY A CHILD.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFRCHMODULE,   ONLY:RECH
C     ------------------------------------------------------------------
C1------ZERO OUT RECHARGE
      RECH(J,I) = 0.
C2------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2RCH
      
C***********************************************************************

C-----VERSION 2.0 15JUNE2013 SGWF2LGR2RES
      SUBROUTINE SGWF2LGR2RES(J,I)
C     ******************************************************************
C     ZERO OUT RESERVOIR LEAKAGE IN THE PARENT AREA THAT IS COVERED BY 
C     A CHILD.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFRESMODULE,   ONLY:CRES
C     ------------------------------------------------------------------
C1------ZERO OUT RESERVOIR LEAKAGE
      CRES(J,I) = 0.
C2------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2RES
      
C***********************************************************************

C-----VERSION 1.0 13AUGUST2010 SGWF2LGR2WETCHK
      SUBROUTINE SGWF2LGR2WETCHK(J,I,K,HTMP,IHFLG)
C     ******************************************************************
C     CALL EITHER A PARENT OR CHILD ROUTINE TO CHECK IF FOR REWETTING
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE LGRMODULE,   ONLY:ISCHILD
C     ------------------------------------------------------------------
C
C1------FIRST CHECK IF IT IS A PARENT OR CHILD GRID AND CALL ROUTINE
C1------IF PARENT ONLY CALL IF REWETTING FROM HORIZONTAL CONNECTIONS      
      IF(ISCHILD .GT. 0) THEN
        CALL SGWF2LGR2WETCHKC(J,I,K,HTMP,IHFLG) 
      ELSE IF(IHFLG.GT.0) THEN
        CALL SGWF2LGR2WETCHKP(J,I,K,HTMP) 
      END IF
C2------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2WETCHK
      
C***********************************************************************

C-----VERSION 1.0 13AUGUST2010 SGWF2LGR2WETCHKC
      SUBROUTINE SGWF2LGR2WETCHKC(J,I,K,HTMP,IHFLG)
C     ******************************************************************
C     DETERMINE THE MAXIMUM GHOST-NODE HEAD ATTACHED TO THIS CELL.
C     VALUE IS PASSED BACK IN HTMP AND USED TO CHECK FOR REWETTING.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE LGRMODULE,   ONLY:IBFLG,NCON,ICBOUND,IFACEGN,GNHEAD
C     ------------------------------------------------------------------
C
C1------FIRST CHECK IF CELL IS AN INTERFACE CELL.  SKIP IF NOT.
      IF(ICBOUND(J,I,K) .EQ. IBFLG ) THEN
C2------FIND LOCATION OF GHOST-NODE CONNECTION BY SUMMING INDEX VALUE
        IB = COUNT(ICBOUND(1:J,1:I,1:K) .EQ. IBFLG)
C3------LOOP THROUGH ALL GHOST-NODE CONNECTIONS AND USE HIGHEST VALUE
C3------BUT SEPARATE CHECKS FOR VERTICAL AND HORIZONTAL CONNECTIONS
        DO NC=1,NCON(IB)
          IF(IHFLG .EQ. 0 .AND. IFACEGN(NC,IB) .EQ. 6) THEN
            HTMP=GNHEAD(NC,IB)
          ELSE IF(IHFLG .EQ. 1 .AND. IFACEGN(NC,IB) .NE. 6) THEN   
            HTMP=MAX(HTMP,GNHEAD(NC,IB))
          END IF
        END DO
      END IF
C4------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2WETCHKC
      
C***********************************************************************

C-----VERSION 1.0 13AUGUST2010 SGWF2LGR2WETCHKP
      SUBROUTINE SGWF2LGR2WETCHKP(J,I,K,HTMP)
C     ******************************************************************
C     DETERMINE THE MAXIMUM CHILD HEAD ATTACHED TO THIS CELL.
C     VALUE IS PASSED BACK IN HTMP AND USED TO CHECK FOR REWETTING.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE LGRMODULE,   ONLY:NGRDS,LGRDAT
C     ------------------------------------------------------------------
C
C1------FIRST CHECK IF CELL IS AN INTERFACE CELL.  SKIP IF NOT.
      GRIDS: DO LG=2,NGRDS    !swm: will need modification for generations
        DO IB = 1,LGRDAT(LG)%NBNODES
          IF(LGRDAT(LG)%KPLC(1,IB) .EQ. K .AND. 
     1       LGRDAT(LG)%IPLC(1,IB) .EQ. I .AND.
     2       LGRDAT(LG)%JPLC(1,IB) .EQ. J) THEN
C3-------FOUND AN INTERFACE CELL.  LOOP THROUGH ALL GHOST-NODE 
C3-------CONNECTIONS AND USE HIGHEST VALUE.  EXIT LOOPS WHEN DONE
            DO NC=1,LGRDAT(LG)%NCON(IB)
              IF (LGRDAT(LG)%IFACEGN(NC,IB) .NE. 6) HTMP=MAX(HTMP,
     1            LGRDAT(LG)%GNHEAD(NC,IB))
            END DO
            EXIT GRIDS
          END IF
        END DO
      END DO GRIDS
C4------RETURN
      RETURN
      END SUBROUTINE SGWF2LGR2WETCHKP

      end module GwfLgrSubsModule
