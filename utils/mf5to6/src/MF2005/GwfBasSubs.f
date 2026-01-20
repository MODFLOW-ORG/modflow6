      module GwfBasSubs
        use ConstantsModule,     only: MAXCHARLEN, LENBIGLINE
        use ConstantsPHMFModule, only: FCDATAIN, FCDATAOUT, FCDATABOUT,
     &                                 FCDATABIN, FCINPUT,
     &                                 HDRYDEFAULT, HNOFLODEFAULT
        use ChdObsWriterModule, only: ChdObsWriterType,
     &                                createChdObsWriter
        use ChdPackageWriterModule, only: ChdPackageWriterType
        use DisWriterModule, only: DisWriterType
        use FileTypeModule, only: FileType
        use global, only: iout
        use GlobalVariablesModule, only: echo
        use GwfBasModule, only: SGWF2BAS7PNT, SGWF2BAS7PSV
        use MessageModule, only: write_message_centered
        use ModelModule, only: ModelType
        use ObsWriterModule, only: ObsWriterType
        use OpenSpecModule, only: ACCESS, ACTION, FORM
        use PackageWriterModule, only : PackageWriterType
        use SimPHMFModule, only: store_error, store_note, store_warning, 
     &                           ustop
        use TdisWriterModule, only: TdisWriterType
        use UtilitiesModule, only: ConstantReal2D, Write1Drel,
     &                             Write2Drel
        use utl7module, only: U1DREL, U2DREL, USTOPx, URDCOM, URWORD,
     &                        upcase, ulaprwc, u2dint
        use WelPackageWriterModule, only: WelPackageWriterType

        integer :: iprnbas
      contains

C###############################################################################
      SUBROUTINE GWF2BAS7AR(INUNIT,CUNIT,PROGNAM,IUDIS,IUZON,IUMLT,
     &              MAXUNIT,IUOC,HEADNG,IUPVAL,MFVNAM,model)
C     ******************************************************************
C     Allocate and Read for GWF Basic Package
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      use GlobalVariablesModule, only: NIUNIT
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,NPER,NBOTM,NCNFBD,ITMUNI,
     1                     LENUNI,IXSEC,ITRSS,INBAS,IFREFM,NODES,
     2                     MXITER,IUNIT,HNEW,LBOTM,LAYCBD,LAYHDT,
     3                     LAYHDS,PERLEN,NSTP,TSMULT,ISSFLG,DELR,DELC,
     4                     BOTM,HOLD,IBOUND,CR,CC,CV,HCOF,RHS,BUFF,STRT,
     5                     DDREF, constantdelr, constantdelc,
     6                     cbcfilename
      USE PARAMMODULE,ONLY:MXPAR,MXCLST,MXINST,ICLSUM,IPSUM,
     1                     INAMLOC,NMLTAR,NZONAR,NPVAL,
     2                     B,IACTIVE,IPLOC,IPCLST,PARNAM,PARTYP,
     3                     ZONNAM,MLTNAM,INAME, AllocateParamScalars
      USE GWFBASMODULE,ONLY:MSUM,IHEDFM,IHEDUN,IDDNFM,IDDNUN,IBOUUN,
     1                      LBHDSV,LBDDSV,LBBOSV,IBUDFL,ICBCFL,IHDDFL,
     2                      IAUXSV,IBDOPT,IPRTIM,IPEROC,ITSOC,ICHFLG,
     3                      IDDREF,IDDREFNEW,DELT,PERTIM,TOTIM,HNOFLO,
     4                      HDRY,STOPER,CHEDFM,CDDNFM,CBOUFM,VBVL,VBNM,
     5                      AllocateGwfBasScalars
C
      implicit none
      integer :: inunit, iudis, iuzon, iumlt, maxunit, iuoc, iupval, i,
     &           ipbuff, iphnew, istart, istop, j, k, kk, lloc, name, n
      real :: zero
      type(ModelType), pointer, intent(inout) :: model
      character(len=*), intent(in) :: PROGNAM
      double precision :: constantval
      logical :: constant
      integer :: iprn, knew
      character(len=500) :: msg
      CHARACTER*4 CUNIT(NIUNIT)
      CHARACTER*80 HEADNG(2)
      CHARACTER*(*) MFVNAM
      CHARACTER*200 LINE
C
      DOUBLE PRECISION :: HNF, r
      CHARACTER*24 ANAME(2)
      DATA ANAME(1) /'          BOUNDARY ARRAY'/
      DATA ANAME(2) /'            INITIAL HEAD'/
C     ------------------------------------------------------------------
C1------Allocate scalar variables, which makes it possible for multiple
C1------grids to be defined.
      ! Allocation of scalar variables is now done by
      ! AllocateGlobalScalars routine when ModelType is initialized.
      MXITER=1
      if (.not. associated(iunit)) ALLOCATE(IUNIT(NIUNIT))
C
      call AllocateParamScalars()
!      ALLOCATE(ICLSUM,IPSUM,INAMLOC,NMLTAR,NZONAR,NPVAL)
      if (.not. associated(b)) ALLOCATE (B(MXPAR))
      if (.not. associated(iactive)) ALLOCATE (IACTIVE(MXPAR))
      if (.not. associated(iploc)) ALLOCATE (IPLOC(4,MXPAR))
      if (.not. associated(ipclst)) ALLOCATE (IPCLST(14,MXCLST))
      if (.not. associated(parnam)) ALLOCATE (PARNAM(MXPAR))
      if (.not. associated(partyp)) ALLOCATE (PARTYP(MXPAR))
      if (.not. associated(iname)) ALLOCATE (INAME(MXINST))
C
      call AllocateGwfBasScalars()
      constantdelr = 0.0
      constantdelc = 0.0
      cbcfilename = ''
      HDRY=hdrydefault
      STOPER=0.0
      !
      write(*,*)'Processing BAS package input...'
C
C2------Open all files in name file.
      CALL SGWF2BAS7OPEN(INUNIT,IUNIT,CUNIT,NIUNIT,
     &                 PROGNAM,INBAS,MAXUNIT,MFVNAM,model)
C
C3------PRINT A MESSAGE IDENTIFYING THE BASIC PACKAGE.
      WRITE(IOUT,1)INBAS
    1 FORMAT(1X,/1X,'BAS -- BASIC PACKAGE, VERSION 7, 5/2/2005',
     2' INPUT READ FROM UNIT ',I4)
C
C3A-----SHOW PRECISION OF VARIABLES
      IPBUFF=PRECISION(BUFF)
      IPHNEW=PRECISION(HNEW)
      WRITE(IOUT,*)
      IF(IPBUFF.NE.IPHNEW) THEN
        WRITE(IOUT,*) 'MF5to6 was compiled using mixed precision'
        WRITE(IOUT,*) 'Precision of REAL variables:',IPBUFF
        WRITE(IOUT,*) 'Precision of DOUBLE PRECISION variables:',IPHNEW
      ELSE
        WRITE(IOUT,*) 'MF5to6 was compiled using uniform precision'
        WRITE(IOUT,*)
     1   'Precision of REAL and DOUBLE PRECISION variables:',IPBUFF
      END IF
C
C4------Initialize parameter definition variables.
      IPSUM=0
      ICLSUM=0
      INAMLOC=1
      DO 10 N=1,MXPAR
        PARNAM(N)=' '
        PARTYP(N)=' '
        IPLOC(1,N)=0
        IPLOC(2,N)=0
        IACTIVE(N)=0
   10 CONTINUE
C
C5------Allocate and read discretization data.
      CALL SGWF2BAS7ARDIS(IUDIS,IOUT,model)
      NODES=NCOL*NROW*NLAY
      !
      ! -- Now that laycbd has been read, set up layer indices
      !    for fully 3-d discretization.
      call model%InitializeLayptr()
C
C6------Allocate space for global arrays except discretization data.
      ALLOCATE (HNEW(NCOL,NROW,NLAY))
      ALLOCATE (HOLD(NCOL,NROW,NLAY))
      ALLOCATE (IBOUND(NCOL,NROW,NLAY))
      ALLOCATE (CR(NCOL,NROW,NLAY))
      ALLOCATE (CC(NCOL,NROW,NLAY))
      ALLOCATE (CV(NCOL,NROW,NLAY))
      ALLOCATE (HCOF(NCOL,NROW,NLAY))
      ALLOCATE (RHS(NCOL,NROW,NLAY))
      ALLOCATE (BUFF(NCOL,NROW,NLAY))
      ALLOCATE (STRT(NCOL,NROW,NLAY))
      DDREF=>STRT
      ALLOCATE (LAYHDT(NLAY))
      ALLOCATE (LAYHDS(NLAY))
C
C7------Initialize head-dependent thickness indicator to code that
C7------indicates layer is undefined.
      DO 100 I=1,NLAY
        LAYHDT(I)=-1
        LAYHDS(I)=-1
  100 CONTINUE
      WRITE(IOUT,'(//)')
C
C8------Read BAS Package file.
C8A-----READ AND PRINT COMMENTS.  SAVE THE FIRST TWO COMMENTS IN HEADNG.
      HEADNG(1)=' '
      HEADNG(2)=' '
      WRITE(IOUT,*)
      READ(INBAS,'(A)') LINE
      IF(LINE(1:1).NE.'#') GO TO 20
      HEADNG(1)=LINE(1:80)
      WRITE(IOUT,'(1X,A)') HEADNG(1)
      READ(INBAS,'(A)') LINE
      IF(LINE(1:1).NE.'#') GO TO 20
      HEADNG(2)=LINE(1:80)
      WRITE(IOUT,'(1X,A)') HEADNG(2)
      CALL URDCOM(INBAS,IOUT,LINE)
C
C8B-----LOOK FOR OPTIONS IN THE FIRST ITEM AFTER THE HEADING.
   20 IXSEC=0
      ICHFLG=0
      IFREFM=0
      IPRTIM=0
      STOPER=0.0
      LLOC=1
   25 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INBAS)
      IF(LINE(ISTART:ISTOP).EQ.'XSECTION') THEN
         IXSEC=1
      ELSE IF(LINE(ISTART:ISTOP).EQ.'CHTOCH') THEN
         ICHFLG=1
      ELSE IF(LINE(ISTART:ISTOP).EQ.'FREE') THEN
         IFREFM=1
         WRITE(IOUT,26)
   26    FORMAT (1X,'THE FREE FORMAT OPTION HAS BEEN SELECTED')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'PRINTTIME') THEN
         IPRTIM=1
         WRITE(IOUT,7)
    7    FORMAT(1X,'THE PRINTTIME OPTION HAS BEEN SELECTED')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'STOPERROR') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,N,STOPER,IOUT,INBAS)
         WRITE(IOUT,8) STOPER
    8    FORMAT(1X,'When solver convergence criteria are not met,',/
     1    1X,'execution will continue unless the budget percent',/
     2    1X,'discrepancy is greater than:',F10.4)
      END IF
      IF(LLOC.LT.200) GO TO 25
C
C8C-----PRINT A MESSAGE SHOWING OPTIONS.
      IF(IXSEC.NE.0) WRITE(IOUT,61)
   61 FORMAT(1X,'CROSS SECTION OPTION IS SPECIFIED')
      IF(ICHFLG.NE.0) WRITE(IOUT,62)
   62 FORMAT(1X,'CALCULATE FLOW BETWEEN ADJACENT CONSTANT-HEAD CELLS')
C
C8D-----INITIALIZE TOTAL ELAPSED TIME COUNTER STORAGE ARRAY COUNTER
C8D-----AND CALCULATE NUMBER OF CELLS.
      TOTIM=0.
C
C8E-----READ BOUNDARY ARRAY(IBOUND).
      IF(IXSEC.EQ.0) THEN
         DO 280 K=1,NLAY
         KK=K
         CALL U2DINT(IBOUND(:,:,KK),ANAME(1),NROW,NCOL,KK,INBAS,IOUT)
  280    CONTINUE
      ELSE
         CALL U2DINT(IBOUND(:,:,1),ANAME(1),NLAY,NCOL,-1,INBAS,IOUT)
      END IF
      !
C
C8F-----READ AND PRINT HEAD VALUE TO BE PRINTED FOR NO-FLOW CELLS.
      IF(IFREFM.EQ.0) THEN
         READ(INBAS,'(F10.0)') HNOFLO
      ELSE
         READ(INBAS,*) HNOFLO
      END IF
      HNOFLO = hnoflodefault
      HNF=HNOFLO
      model%NpfWriter%Hnoflo = HNOFLO
      WRITE(IOUT,3) HNOFLO
    3 FORMAT(1X,/1X,'AQUIFER HEAD WILL BE SET TO ',1PG12.5,
     1       ' AT ALL NO-FLOW NODES (IBOUND=0).')
C
C8G-----READ INITIAL HEADS.
      if (.not. model%ConversionDone) then
!       Write BEGIN line for IcData block
        write(model%IcWriter%fileobj%IUnit,285)
  285   format('BEGIN GRIDDATA',/,2x,'STRT LAYERED')
      endif
      IF(IXSEC.EQ.0) THEN
         knew = 0
         DO 300 K=1,NLAY
         KK=K
         knew = knew + 1
         CALL U2DREL(STRT(:,:,KK),ANAME(2),NROW,NCOL,KK,INBAS,IOUT,
     &               constantval,constant,iprn, model%Mf2005Files)
         call ConstantReal2D(ncol, nrow, strt(:,:,kk), constant,
     &                       constantval)
         if (echo) then
           iprnbas = abs(iprn)
         else
           iprnbas = -abs(iprn)
         endif
         ! Write STRT to IC8 file
         write(msg,255)'Initial head for layer ',knew
  255    format(a,i0)
         if (.not. model%ConversionDone) then
           call Write2Drel(model%IcWriter%fileobj%IUnit, nrow, ncol,
     &       STRT(:,:,KK), constant, constantval, msg, .false., iprnbas)
         endif
         ! If layer has an underlying quasi-3d confining bed, write
         ! starting head for it here. Assume it equals head in the
         ! layer. Call store_note with message to this effect.
         if (LAYCBD(k) /= 0) then
           knew = knew + 1
           if (.not. model%ConversionDone) then
             call Write2Drel(model%IcWriter%fileobj%IUnit, nrow, ncol,
     &       STRT(:,:,KK), constant, constantval, msg, .false., iprnbas)
           endif
           write(msg,260)knew,k
           call store_note(msg)
         endif
  260    format('Starting head for MODFLOW 6 model layer ',i0,
     &   ', which replaces a quasi-3D unit, has been assigned equal'
     &   ' to starting head for layer ',i0,' of original model.')
  300    CONTINUE
      ELSE
         ! for cross-section model
         CALL U2DREL(STRT(:,:,1),ANAME(2),NLAY,NCOL,-1,INBAS,IOUT,
     &               constantval,constant,iprn,model%Mf2005Files)
         call ConstantReal2D(ncol, nrow, strt(:,:,1), constant,
     &                       constantval)
         if (echo) then
           iprnbas = abs(iprn)
         else
           iprnbas = -abs(iprn)
         endif
         ! Write STRT to IC8 file
         msg = 'Initial head for cross section'
         if (.not. model%ConversionDone) then
           do k=1,nlay
             call Write2Drel(model%IcWriter%fileobj%IUnit, nrow, ncol,
     &       STRT(:,:,k), constant, constantval, msg, .false., iprnbas)
           enddo
         endif
      END IF
      if (.not. model%ConversionDone) then
!       Write END line for IcData block
        write(model%IcWriter%fileobj%IUnit,305)
  305   format('END GRIDDATA')
      endif

      ! Now that STRT has been read, have IbChdWriter extract CH cells from IBOUND
      model%IbChdWriter%ibnd => IBOUND
      call model%IbChdWriter%ProcessIbound()
      if (model%IbChdWriter%IbChdList%Count() == 0) then
        model%IbChdWriter%ObsActive = .false.
        model%IbChdWriter%Active = .false.
      endif
C
C9------COPY INITIAL HEADS FROM STRT TO HNEW.
      DO 402 K=1,NLAY
        DO 401 I=1,NROW
          DO 400 J=1,NCOL
            HNEW(J,I,K)=STRT(J,I,K)
            IF(IBOUND(J,I,K).EQ.0) HNEW(J,I,K)=HNF
  400     CONTINUE
  401   continue
  402 continue
C
C10-----SET UP OUTPUT CONTROL.
      CALL SGWF2BAS7I(NLAY,IUNIT(IUOC),IOUT,IFREFM,NIUNIT,model)
C
C11-----INITIALIZE VOLUMETRIC BUDGET ACCUMULATORS TO ZERO.
  590 ZERO=0.
      DO 600 I=1,NIUNIT
      DO 600 J=1,4
      VBVL(J,I)=ZERO
  600 CONTINUE
C
C12-----Allocate and read Zone and Multiplier arrays
      CALL SGWF2BAS7ARMZ(IUNIT(IUZON),IUNIT(IUMLT))
C
C13-----READ PARAMETER VALUES FILE.
      CALL SGWF2BAS7ARPVAL(IUPVAL)
C
C14-----SAVE POINTERS TO DATA AND RETURN.
      CALL SGWF2BAS7PSV(model%IGrid)
      RETURN
      END subroutine GWF2BAS7AR

C###############################################################################

      SUBROUTINE SGWF2BAS7ARDIS(IUDIS,IOUT,model)
C     *****************************************************************
C     ALLOCATE AND READ DIS DATA
C     *****************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,NPER,NBOTM,NCNFBD,ITMUNI,
     1                     LENUNI,IUNIT,LBOTM,LAYCBD,ITRSS,
     3                     PERLEN,NSTP,TSMULT,ISSFLG,DELR,DELC,BOTM
C
      type(ModelType), intent(inout) :: model
      CHARACTER*200 LINE
      CHARACTER*24 ANAME(5)
      DATA ANAME(1) /'                    DELR'/
      DATA ANAME(2) /'                    DELC'/
      DATA ANAME(3) /'TOP ELEVATION OF LAYER 1'/
      DATA ANAME(4) /'  MODEL LAYER BOTTOM EL.'/
      DATA ANAME(5) /'BOT. EL. OF QUASI-3D BED'/
      ! local variables
      logical :: constant
      character(len=500) :: msg
      double precision :: r, constantval
  910 format(a,1x,i0)
C     ------------------------------------------------------------------
C

C1------Check for existence of discretization file
      INDIS=IUNIT(IUDIS)
      IF(INDIS.LE.0) THEN
         WRITE(IOUT,*) ' DIS file must be specified for MODFLOW to run'
         CALL USTOPx(' ')
      END IF
      WRITE(IOUT,11) INDIS
   11 FORMAT(1X,/1X,'DISCRETIZATION INPUT DATA READ FROM UNIT ',I4)
C
C
C2------Read comments and the first line following the comments.
      CALL URDCOM(INDIS,IOUT,LINE)
C
C3------Get the number of layers, rows, columns, stress periods,
C3------ITMUNI, and LENUNI from the line.
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NLAY,R,IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NROW,R,IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NCOL,R,IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPER,R,IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITMUNI,R,IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,LENUNI,R,IOUT,INDIS)
C
C4------PRINT # OF LAYERS, ROWS, COLUMNS AND STRESS PERIODS.
      WRITE(IOUT,15) NLAY,NROW,NCOL
   15 FORMAT(1X,I4,' LAYERS',I10,' ROWS',I10,' COLUMNS')
      model%nrow = NROW
      model%ncol = NCOL
      WRITE(IOUT,20) NPER
   20 FORMAT(1X,I4,' STRESS PERIOD(S) IN SIMULATION')
      call model%TdisWriter%Alloc(nper)
C
C5------SELECT AND PRINT A MESSAGE SHOWING TIME UNIT.
      IF(ITMUNI.LT.0 .OR. ITMUNI.GT.5) ITMUNI=0
      IF(ITMUNI.EQ.0) THEN
         WRITE(IOUT,30)
   30    FORMAT(1X,'MODEL TIME UNIT IS UNDEFINED')
         call model%AssignTimeUnit('UNDEFINED')
      ELSE IF(ITMUNI.EQ.1) THEN
         WRITE(IOUT,40)
   40    FORMAT(1X,'MODEL TIME UNIT IS SECONDS')
         call model%AssignTimeUnit('SECONDS')
      ELSE IF(ITMUNI.EQ.2) THEN
         WRITE(IOUT,50)
   50    FORMAT(1X,'MODEL TIME UNIT IS MINUTES')
         call model%AssignTimeUnit('MINUTES')
      ELSE IF(ITMUNI.EQ.3) THEN
         WRITE(IOUT,60)
   60    FORMAT(1X,'MODEL TIME UNIT IS HOURS')
         call model%AssignTimeUnit('HOURS')
      ELSE IF(ITMUNI.EQ.4) THEN
         WRITE(IOUT,70)
   70    FORMAT(1X,'MODEL TIME UNIT IS DAYS')
         call model%AssignTimeUnit('DAYS')
      ELSE
         WRITE(IOUT,80)
   80    FORMAT(1X,'MODEL TIME UNIT IS YEARS')
         call model%AssignTimeUnit('YEARS')
      END IF
C
C6------SELECT AND PRINT A MESSAGE SHOWING LENGTH UNIT.
      IF(LENUNI.LT.0 .OR. LENUNI.GT.3) LENUNI=0
      IF(LENUNI.EQ.0) THEN
         WRITE(IOUT,90)
   90    FORMAT(1X,'MODEL LENGTH UNIT IS UNDEFINED')
         call model%AssignLengthUnit('UNDEFINED')
      ELSE IF(LENUNI.EQ.1) THEN
         WRITE(IOUT,91)
   91    FORMAT(1X,'MODEL LENGTH UNIT IS FEET')
         call model%AssignLengthUnit('FEET')
      ELSE IF(LENUNI.EQ.2) THEN
         WRITE(IOUT,93)
   93    FORMAT(1X,'MODEL LENGTH UNIT IS METERS')
         call model%AssignLengthUnit('METERS')
      ELSE IF(LENUNI.EQ.3) THEN
         WRITE(IOUT,95)
   95    FORMAT(1X,'MODEL LENGTH UNIT IS CENTIMETERS')
         call model%AssignLengthUnit('CENTIMETERS')
      END IF
C
C7------ALLOCATE LAYER FLAGS.
      ALLOCATE(LBOTM(NLAY))
      ALLOCATE(LAYCBD(NLAY))
C
C8------Read confining bed information
      READ(INDIS,*) (LAYCBD(K),K=1,NLAY)
      LAYCBD(NLAY)=0
      WRITE(IOUT,*) ' Confining bed flag for each layer:'
      WRITE(IOUT,'(20I4)') (LAYCBD(K),K=1,NLAY)
C
C9------Count confining beds, setup the pointer to each layer's
C9------bottom array (LBOTM), and setup LAYCBD to be the confining
C9------bed number for each layer.
      NCNFBD=0
      DO 100 K=1,NLAY
      LBOTM(K)=K+NCNFBD
      IF(LAYCBD(K).NE.0) THEN
         NCNFBD=NCNFBD+1
         LAYCBD(K)=NCNFBD
      END IF
  100 CONTINUE
      if (NCNFBD > 0) then
        msg = 'MODFLOW 6 does not support quasi-3D confining' //
     &        ' units. MF5to6 generates an active layer to' //
     &        ' replace each quasi-3D confining unit, and default' //
     &        ' values are provided for hydraulic properties' //
     &        ' of the new layers. These values may need' //
     &        ' to be adjusted by the modeler to approximate the' //
     &        ' original model.'
        call store_note(msg)
      endif
      NBOTM=NLAY+NCNFBD
      model%nbotm = NBOTM
C
C10-----Allocate space for discretization arrays
C10-----Note that NBOTM+1 arrays are allocated for BOTM
C10-----because BOTM(J,I,0) contains the top elevation of layer 1.
      ALLOCATE (DELR(NCOL))
      ALLOCATE (DELC(NROW))
      ALLOCATE (BOTM(NCOL,NROW,0:NBOTM))
      ALLOCATE (PERLEN(NPER),NSTP(NPER),TSMULT(NPER),ISSFLG(NPER))
C
C11-----Read the DELR and DELC arrays.
      CALL U1DREL(DELR,ANAME(1),NCOL,INDIS,IOUT, constantval,
     &    constant, iprn)
      !
      CALL U1DREL(DELC,ANAME(2),NROW,INDIS,IOUT,constantval,
     &    constant, iprn)
C
C12-----Read the top elevation of layer 1.
      CALL U2DREL(BOTM(:,:,0),ANAME(3),NROW,NCOL,0,INDIS,IOUT,
     &    constantval, constant, iprn)
C
C13-----Read the bottom elevations.
      DO 120 K=1,NLAY
        KK=K
        CALL U2DREL(BOTM(:,:,LBOTM(K)),ANAME(4),NROW,NCOL,KK,INDIS,
     &      IOUT, constantval, constant, iprn)
        IF(LAYCBD(K).NE.0) then
          CALL U2DREL(BOTM(:,:,LBOTM(K)+1),ANAME(5),
     1          NROW,NCOL,KK,INDIS,IOUT, constantval, constant, iprn)
          ! Write array to DIS file
          if (echo) then
            iprnbas = abs(iprn)
          else
            iprnbas = -abs(iprn)
          endif
        endif
  120 CONTINUE
C
C14-----READ AND WRITE LENGTH OF STRESS PERIOD, NUMBER OF TIME STEPS,
C14-----TIME STEP MULTIPLIER, AND STEADY-STATE FLAG..
      WRITE(IOUT,161)
  161 FORMAT(1X,//1X,'STRESS PERIOD     LENGTH       TIME STEPS',
     1            '     MULTIPLIER FOR DELT    SS FLAG',/1X,76('-'))
      ISS=0
      ITR=0
      DO 200 N=1,NPER
      READ(INDIS,'(A)') LINE
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,PERLEN(N),IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NSTP(N),R,IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,TSMULT(N),IOUT,INDIS)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,INDIS)
      model%TdisWriter%StressPeriods(n)%nstp = nstp(n)
      model%TdisWriter%StressPeriods(n)%perlen = PERLEN(n)
      model%TdisWriter%StressPeriods(n)%tsmult = TSMULT(n)
      model%TdisWriter%StressPeriods(n)%sstr = LINE(ISTART:ISTOP)
      IF (LINE(ISTART:ISTOP).EQ.'TR') THEN
         ISSFLG(N)=0
         ITR=1
         if (.not. model%StoWriter%Active) then
            model%StoWriter%Active = .true.
            call model%StoWriter%AllocateArrays()
            call model%StoWriter%InitializeFile(
     &        model%StoWriter%fileobj%FName, 'STO6')
            call model%Mf6Files%AddFile(
     &        model%StoWriter%fileobj%FName,
     &        'STO6', model%StoWriter%fileobj%IUnit,
     &        model%StoWriter%fileobj%FCode,
     &        model%StoWriter%fileobj%PkgName)
         endif
         model%StoWriter%Transient(n) = .true.
      ELSE IF (LINE(ISTART:ISTOP).EQ.'SS') THEN
         ISSFLG(N)=1
         ISS=1
      ELSE
         WRITE(IOUT,162)
  162    FORMAT(' SSFLAG MUST BE EITHER "SS" OR "TR"',
     1      ' -- STOP EXECUTION (SGWF2BAS7ARDIS)')
         CALL USTOPx(' ')
      END IF
      WRITE (IOUT,163) N,PERLEN(N),NSTP(N),TSMULT(N),LINE(ISTART:ISTOP)
  163 FORMAT(1X,I8,1PG21.7,I7,0PF25.3,A11)
C
C15-----STOP IF NSTP LE 0, PERLEN EQ 0 FOR TRANSIENT STRESS PERIODS,
C15-----TSMULT LE 0, OR PERLEN LT 0..
      IF(NSTP(N).LE.0) THEN
         WRITE(IOUT,164)
  164    FORMAT(1X,/1X,
     1  'THERE MUST BE AT LEAST ONE TIME STEP IN EVERY STRESS PERIOD')
         CALL USTOPx(' ')
      END IF
      ZERO=0.
      IF(PERLEN(N).EQ.ZERO .AND. ISSFLG(N).EQ.0) THEN
         WRITE(IOUT,165)
  165    FORMAT(1X,/1X,
     1  'PERLEN MUST NOT BE 0.0 FOR TRANSIENT STRESS PERIODS')
         CALL USTOPx(' ')
      END IF
      IF(TSMULT(N).LE.ZERO) THEN
         WRITE(IOUT,170)
  170    FORMAT(1X,/1X,'TSMULT MUST BE GREATER THAN 0.0')
         CALL USTOPx(' ')
      END IF
      IF(PERLEN(N).LT.ZERO) THEN
         WRITE(IOUT,175)
  175    FORMAT(1X,/1X,
     1  'PERLEN CANNOT BE LESS THAN 0.0 FOR ANY STRESS PERIOD')
         CALL USTOPx(' ')
      END IF
  200 CONTINUE
C
C16-----Assign ITRSS.
      IF(ISS.EQ.0 .AND. ITR.NE.0) THEN
         ITRSS=1
         WRITE(IOUT,270)
  270    FORMAT(/,1X,'TRANSIENT SIMULATION')
      ELSE IF(ISS.NE.0 .AND. ITR.EQ.0) THEN
         ITRSS=0
         WRITE(IOUT,275)
  275    FORMAT(/,1X,'STEADY-STATE SIMULATION')
      ELSE
         ITRSS=-1
         WRITE(IOUT,280)
  280    FORMAT(/,1X,'COMBINED STEADY-STATE AND TRANSIENT SIMULATION')
      END IF
C
C17-----RETURN.
      RETURN
      END SUBROUTINE SGWF2BAS7ARDIS

C###############################################################################

      SUBROUTINE SGWF2BAS7I(NLAY,INOC,IOUT,IFREFM,NIUNIT,model)
C     ******************************************************************
C     SET UP OUTPUT CONTROL.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFBASMODULE, ONLY: IHEDFM,IDDNFM,IHEDUN,IDDNUN,IPEROC,ITSOC,
     1                        CHEDFM,CDDNFM,IBDOPT,LBHDSV,LBDDSV,
     2                        IBOUUN,LBBOSV,CBOUFM,IAUXSV,IOFLG,
     3                        VBVL,VBNM,IDDREF,IDDREFNEW
      type(ModelType), intent(inout) :: model
      CHARACTER*200 LINE
      character(len=LENBIGLINE) :: msg
      double precision :: r
      type(FileType), pointer :: filtyp
C     ------------------------------------------------------------------
C
C1-----ALLOCATE SPACE FOR IOFLG, VBVL, AND VBNM ARRAYS.
      ALLOCATE (IOFLG(NLAY,5))
      ALLOCATE (VBVL(4,NIUNIT))
      ALLOCATE (VBNM(NIUNIT))
      IDDREF=0
      IDDREFNEW=0
C
C1------ASSIGN DEFAULT VALUES.
      CHEDFM=' '
      CDDNFM=' '
      CBOUFM='(20I4)'
      IHEDFM=0
      IDDNFM=0
      IHEDUN=0
      IDDNUN=0
      IBOUUN=0
      IBDOPT=1
      LBHDSV=0
      LBDDSV=0
      LBBOSV=0
      IAUXSV=0
C
C2------TEST OUTPUT CONTROL INPUT UNIT TO SEE IF OUTPUT CONTROL IS
C2------ACTIVE.
      IF(INOC.LE.0) THEN
C
C2A-----OUTPUT CONTROL IS INACTIVE. PRINT A MESSAGE LISTING DEFAULTS.
         WRITE(IOUT, 41)
   41    FORMAT(1X,/1X,'DEFAULT OUTPUT CONTROL',/1X,
     1   'THE FOLLOWING OUTPUT COMES AT THE END OF EACH STRESS PERIOD:')
         WRITE(IOUT, 42)
   42    FORMAT(1X,'TOTAL VOLUMETRIC BUDGET')
         WRITE(IOUT, 43)
   43    FORMAT(1X,10X,'HEAD')
C
C2B-----SET DEFAULT FLAGS IN IOFLG SO THAT HEAD IS PRINTED FOR
C2B-----EVERY LAYER.
         DO 80 K=1,NLAY
         IOFLG(K,1)=1
         IOFLG(K,2)=0
         IOFLG(K,3)=0
         IOFLG(K,4)=0
         IOFLG(K,5)=0
   80    CONTINUE
         GO TO 1000
      END IF
C
      model%OcWriter%Active = .true.
      model%OcWriter%ModelBasename = model%BaseName
      call model%InitializeMF6File('OC6', FCINPUT)
C3------OUTPUT CONTROL IS ACTIVE.  READ FIRST RECORD AND DECODE FIRST
C3------WORD.  MUST USE URWORD IN CASE FIRST WORD IS ALPHABETIC.
      CALL URDCOM(INOC,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
C
C4------TEST FOR NUMERIC OUTPUT CONTROL.  FIRST WORD WILL NOT BE
C4------"PERIOD", "HEAD", "DRAWDOWN", OR "COMPACT".
      IF(LINE(ISTART:ISTOP).NE.'PERIOD' .AND. LINE(ISTART:ISTOP).NE.
     1     'HEAD' .AND. LINE(ISTART:ISTOP).NE.'DRAWDOWN' .AND.
     2     LINE(ISTART:ISTOP).NE.'COMPACT' .AND.
     3     LINE(ISTART:ISTOP).NE.'IBOUND') THEN
C4A-----NUMERIC OUTPUT CONTROL.  DECODE THE INITIAL RECORD ACCORDINGLY.
         WRITE(IOUT,102)
  102    FORMAT(1X,/1X,'OUTPUT CONTROL IS SPECIFIED EVERY TIME STEP')
         IF(IFREFM.EQ.0) THEN
            READ(LINE,'(4I10)') IHEDFM,IDDNFM,IHEDUN,IDDNUN
         ELSE
            LLOC=1
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDFM,R,IOUT,INOC)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNFM,R,IOUT,INOC)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDUN,R,IOUT,INOC)
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNUN,R,IOUT,INOC)
         END IF
         WRITE(IOUT,103) IHEDFM,IDDNFM
  103    FORMAT(1X,'HEAD PRINT FORMAT CODE IS',I4,
     1     '    DRAWDOWN PRINT FORMAT CODE IS',I4)
         WRITE(IOUT,104) IHEDUN,IDDNUN
  104    FORMAT(1X,'HEADS WILL BE SAVED ON UNIT ',I4,
     1     '    DRAWDOWNS WILL BE SAVED ON UNIT ',I4)
         IPEROC=-1
         ITSOC=-1
         model%OcWriter%IHEDFM = IHEDFM
         model%OcWriter%IDDNFM = IDDNFM
         model%OcWriter%IHEDUN = IHEDUN
         model%OcWriter%IDDNUN = IDDNUN
         !
         ! Find file associated with IHEDUN
         if (IHEDUN > 0) then
           filtyp => model%Mf6Files%GetFileByUnit(IHEDUN)
           model%OcWriter%HdSvFil = filtyp%FName
         endif
         !
         ! Find file associated with IDDNUN
         if (IDDNUN > 0) then
           filtyp => model%Mf6Files%GetFileByUnit(IDDNUN)
           model%OcWriter%HdSvFil = filtyp%FName
         endif
      ELSE
C4B-----ALPHABETIC OUTPUT CONTROL.  CALL MODULE TO READ INITIAL RECORDS.
         CALL SGWF2BAS7J(INOC,IOUT,LINE,LLOC,ISTART,ISTOP,model)
      END IF

      ! -- Store IHEDFM and IDDNFM in OutputControlWriter
      model%OcWriter%IHEDFM = IHEDFM
      model%OcWriter%IDDNFM = IDDNFM
C
C5------RETURN.
 1000 RETURN
      END SUBROUTINE SGWF2BAS7I

C###############################################################################

      SUBROUTINE SGWF2BAS7J(INOC,IOUT,LINE,LLOC,ISTART,ISTOP,model)
C     ******************************************************************
C     READ INITIAL ALPHABETIC OUTPUT CONTROL RECORDS.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GWFBASMODULE, ONLY: IHEDFM,IDDNFM,IHEDUN,IDDNUN,IPEROC,ITSOC,
     1                        CHEDFM,CDDNFM,IBDOPT,LBHDSV,LBDDSV,
     2                        IBOUUN,LBBOSV,CBOUFM,IAUXSV,IDDREFNEW
      type(ModelType), intent(inout) :: model
      type(FileType), pointer :: filtyp
C
      CHARACTER*200 LINE, msg
      double precision :: r
C     ------------------------------------------------------------------
C
C1------ALPHABETIC OUTPUT CONTROL.  WRITE MESSAGE AND SET INITIAL VALUES
C1------FOR IPEROC AND ITSOC.
      WRITE(IOUT,91)
   91 FORMAT(1X,/1X,'OUTPUT CONTROL IS SPECIFIED ONLY AT TIME STEPS',
     1    ' FOR WHICH OUTPUT IS DESIRED')
      IPEROC=9999
      ITSOC=9999
C
C2------LOOK FOR ALPHABETIC WORDS:
C2A-----LOOK FOR "PERIOD", WHICH INDICATES THE END OF INITIAL OUTPUT
C2A-----CONTROL DATA.  IF FOUND, DECODE THE PERIOD NUMBER AND TIME
C2A-----STEP NUMBER FOR LATER USE.
  100 IF(LINE(ISTART:ISTOP).EQ.'PERIOD') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPEROC,R,IOUT,INOC)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
         IF(LINE(ISTART:ISTOP).NE.'STEP') GO TO 2000
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ITSOC,R,IOUT,INOC)
         WRITE(IOUT,101) IHEDFM,IDDNFM
  101    FORMAT(1X,'HEAD PRINT FORMAT CODE IS',I4,
     1        '    DRAWDOWN PRINT FORMAT CODE IS',I4)
         WRITE(IOUT,102) IHEDUN,IDDNUN
  102    FORMAT(1X,'HEADS WILL BE SAVED ON UNIT ',I4,
     1        '    DRAWDOWNS WILL BE SAVED ON UNIT ',I4)
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
         IF(LINE(ISTART:ISTOP).EQ.'DDREFERENCE') THEN
           IDDREFNEW=1
         ELSE
           IDDREFNEW=0
         END IF
         GO TO 1000
C
C2B-----LOOK FOR "HEAD PRINT ..." AND "HEAD SAVE ...".  IF
C2B-----FOUND, SET APPROPRIATE FLAGS.
      ELSE IF(LINE(ISTART:ISTOP).EQ.'HEAD') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
         IF(LINE(ISTART:ISTOP).EQ.'PRINT') THEN
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
            IF(LINE(ISTART:ISTOP).NE.'FORMAT') GO TO 2000
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDFM,R,IOUT,INOC)
            model%OcWriter%IHEDFM = IHEDFM
         ELSE IF(LINE(ISTART:ISTOP).EQ.'SAVE') THEN
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
            IF(LINE(ISTART:ISTOP).EQ.'UNIT') THEN
               CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IHEDUN,R,IOUT,
     1            INOC)
               model%OcWriter%IHEDUN = IHEDUN
            ELSE IF(LINE(ISTART:ISTOP).EQ.'FORMAT') THEN
              msg = 'HEAD SAVE FORMAT is not supported in MODFLOW 6.'
              call store_note(msg)
*               CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,INOC)
*               CHEDFM=LINE(ISTART:ISTOP)
*               model%OcWriter%HdSvFmt = CHEDFM
*               WRITE(IOUT,103) CHEDFM
*  103          FORMAT(1X,'HEADS WILL BE SAVED WITH FORMAT: ',A)
*               CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
*               IF(LINE(ISTART:ISTOP).EQ.'LABEL') THEN
*                  LBHDSV=1
*                  model%OcWriter%HdLbl = .true.
*                  WRITE(IOUT,104)
*  104             FORMAT(1X,'SAVED HEADS WILL BE LABELED')
*               END IF
            ELSE
               GO TO 2000
            END IF
         ELSE
            GO TO 2000
         END IF
C
C2C-----LOOK FOR "DRAWDOWN PRINT ..." AND "DRAWDOWN SAVE ...".
C2C-----IF FOUND, SET APPROPRIATE FLAGS
      ELSE IF(LINE(ISTART:ISTOP).EQ.'DRAWDOWN') THEN
        msg = 'Calculation of drawdown is not supported in MODFLOW 6.'
        call store_note(msg)
*         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
*         IF(LINE(ISTART:ISTOP).EQ.'PRINT') THEN
*            CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
*            IF(LINE(ISTART:ISTOP).NE.'FORMAT') GO TO 2000
*            CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNFM,R,IOUT,INOC)
*            model%OcWriter%IDDNFM = IDDNFM
*         ELSE IF(LINE(ISTART:ISTOP).EQ.'SAVE') THEN
*            CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
*            IF(LINE(ISTART:ISTOP).EQ.'UNIT') THEN
*               CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IDDNUN,R,IOUT,
*     1                   INOC)
*               model%OcWriter%IDDNUN = IDDNUN
*            ELSE IF(LINE(ISTART:ISTOP).EQ.'FORMAT') THEN
*               CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,INOC)
*               CDDNFM=LINE(ISTART:ISTOP)
*               model%OcWriter%DdSvFmt = CDDNFM
*               WRITE(IOUT,105) CDDNFM
*  105          FORMAT(1X,'DRAWDOWN WILL BE SAVED WITH FORMAT: ',A)
*               CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
*               IF(LINE(ISTART:ISTOP).EQ.'LABEL') THEN
*                  LBDDSV=1
*                  model%OcWriter%DdLbl = .true.
*                  WRITE(IOUT,106)
*  106             FORMAT(1X,'SAVED DRAWDOWN WILL BE LABELED')
*               END IF
*            ELSE
*               GO TO 2000
*            END IF
*         ELSE
*            GO TO 2000
*         END IF
C
C2D-----LOOK FOR "COMPACT BUDGET FILES" -- "COMPACT" IS SUFFICIENT.
C2D-----IF FOUND, SET APPROPRIATE FLAG.
      ELSE IF(LINE(ISTART:ISTOP).EQ.'COMPACT') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
         IF(LINE(ISTART:ISTOP).EQ.'BUDGET') THEN
            IBDOPT=2
            WRITE(IOUT,107)
  107       FORMAT(1X,
     1      'COMPACT CELL-BY-CELL BUDGET FILES WILL BE WRITTEN')
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
            IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1         LINE(ISTART:ISTOP).EQ.'AUX') THEN
               IAUXSV=1
               WRITE(IOUT,108)
  108          FORMAT(1X,
     1     'AUXILIARY DATA WILL BE SAVED IN CELL-BY-CELL BUDGET FILES')
            END IF
         ELSE
            GO TO 2000
         END IF
C
C2E-----LOOK FOR  "IBOUND SAVE ...".  IF FOUND, SET APPROPRIATE FLAGS.
      ELSE IF(LINE(ISTART:ISTOP).EQ.'IBOUND') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
         IF(LINE(ISTART:ISTOP).EQ.'SAVE') THEN
            CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
            IF(LINE(ISTART:ISTOP).EQ.'UNIT') THEN
               CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IBOUUN,R,IOUT,
     1            INOC)
               model%OcWriter%IBOUUN = IBOUUN
               WRITE(IOUT,111) IBOUUN
  111          FORMAT(1X,'IBOUND WILL BE SAVED ON UNIT ',I4)
            ELSE IF(LINE(ISTART:ISTOP).EQ.'FORMAT') THEN
               CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,INOC)
               CBOUFM=LINE(ISTART:ISTOP)
               model%OcWriter%IbSvFmt = CBOUFM
               WRITE(IOUT,112) CBOUFM
  112          FORMAT(1X,'IBOUND WILL BE SAVED WITH FORMAT: ',A)
               CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
               IF(LINE(ISTART:ISTOP).EQ.'LABEL') THEN
                  LBBOSV=1
                  model%OcWriter%IbLbl = .true.
                  WRITE(IOUT,109)
  109             FORMAT(1X,'SAVED IBOUND WILL BE LABELED')
               END IF
            ELSE
               GO TO 2000
            END IF
         ELSE
            GO TO 2000
         END IF
C
C2F-----ERROR IF UNRECOGNIZED WORD.
      ELSE
         GO TO 2000
      END IF
C
C3------FINISHED READING A RECORD.  READ NEXT RECORD, IGNORING BLANK
C3------LINES.  GO BACK AND DECODE IT.
  110 READ(INOC,'(A)',END=1000) LINE
      IF(LINE.EQ.' ') GO TO 110
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INOC)
      GO TO 100
C
C4------RETURN.
 1000 RETURN
C
C5------ERROR DECODING INPUT DATA.
 2000 WRITE(IOUT,2001) LINE
 2001 FORMAT(1X,/1X,'ERROR READING OUTPUT CONTROL INPUT DATA:'/1X,A80)
      CALL USTOPx(' ')
      END SUBROUTINE SGWF2BAS7J

C###############################################################################

      SUBROUTINE SGWF2BAS7OPEN(INUNIT,IUNIT,CUNIT,NIUNIT,
     1              PROGNAM,INBAS,MAXUNIT,MFVNAM,model)
C     ******************************************************************
C     OPEN FILES.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      type(ModelType), pointer, intent(inout) :: model
      character(len=*) :: PROGNAM
      DIMENSION IUNIT(NIUNIT)
      CHARACTER*4 CUNIT(NIUNIT)
      CHARACTER*7 FILSTAT
      CHARACTER*20 FILACT, FMTARG, ACCARG
      CHARACTER*(*) MFVNAM
      CHARACTER*40 SPACES
      CHARACTER*300 LINE
      character(len=MAXCHARLEN) :: fname
      CHARACTER*20 FILTYP !, filtyp15
      LOGICAL LOP
      character(len=100)  :: msg
      integer :: fcodelocal !, iu15
      class(PackageWriterType), pointer :: PkgWriter => null()
      type(ChdObsWriterType), pointer :: newChdObsWriter => null()
      type(ChdPackageWriterType), pointer :: chdpw => null()
      double precision :: r
      integer :: i, iuchob, npkgs
      !
      ! -- formats
1     format(/,'MODFLOW-2005 name file: ',a)
2     format('MODFLOW 6 base-name: ',a)
C     ---------------------------------------------------------------
C
C1------INITIALIZE CONSTANTS.
      INBAS=0
      NFILE=0
      IOUT=0
      iuchob = 0
      DO 5 I=1,NIUNIT
      IUNIT(I)=0
5     CONTINUE
      SPACES=' '
      LENVER=LEN_TRIM(PROGNAM)
      INDENT=40-(LENVER+8)/2
C
C2------READ A LINE; IGNORE BLANK LINES AND PRINT COMMENT LINES.
10    READ(INUNIT,'(A)',END=1000) LINE
      IF(LINE.EQ.' ') GO TO 10
      IF(LINE(1:1).EQ.'#') THEN
        IF(NFILE.NE.0 .AND. IOUT.NE.0) WRITE(IOUT,'(A)') LINE
        GO TO 10
      END IF
C
C3------DECODE THE FILE TYPE, UNIT NUMBER, AND NAME.
      LLOC=1
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUT,INUNIT)
      FILTYP=LINE(ITYP1:ITYP2)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IU,R,IOUT,INUNIT)
      CALL URWORD(LINE,LLOC,INAM1,INAM2,0,N,R,IOUT,INUNIT)
      IFLEN=INAM2-INAM1+1
      if (iout<1) then
        fname = trim(model%BaseName) // '_conversion_messages.txt'
        iflen = len_trim(fname)
      else
        FNAME = LINE(INAM1:INAM2)
      endif

      INQUIRE(UNIT=IU,OPENED=LOP)
      IF(LOP) THEN
         IF(IOUT.EQ.0) THEN
            WRITE(*,11) FNAME(1:IFLEN),IU
   11       FORMAT(1X,/1X,'CANNOT OPEN ',A,' ON UNIT',I4,
     1              ' BECAUSE UNIT IS ALREADY BEING USED')
         ELSE
            WRITE(IOUT,11) FNAME(1:IFLEN),IU
         END IF
         CALL USTOPx(' ')
      END IF
C
C4------KEEP TRACK OF LARGEST UNIT NUMBER
      IF (IU.GT.MAXUNIT) MAXUNIT = IU
C
C5------SET DEFAULT FILE ATTRIBUTES.
      FMTARG='FORMATTED'
      ACCARG='SEQUENTIAL'
      FILSTAT='UNKNOWN'
      FILACT=' '
C
C6------SPECIAL CHECK FOR 1ST FILE.
      IF(NFILE.EQ.0) THEN
        IF(FILTYP.EQ.'LIST') THEN
          IOUT=IU
          model%iulist = iu
          OPEN(UNIT=IU,FILE=FNAME(1:IFLEN),STATUS='REPLACE',
     1          FORM='FORMATTED',ACCESS='SEQUENTIAL')
          call write_message_centered(PROGNAM, 80, iunit=iout)
          msg = 'Conversion Report'
          call write_message_centered(msg, 80, iunit=iout)
          write(iout,1)trim(model%NameFile2005)
          write(iout,2)trim(model%BaseName)
          if (.not. model%ConversionDone) then
            call model%InitializeMF6Files()
          endif
        ELSE
          WRITE(*,*)
     1       ' FIRST ENTRY IN NAME FILE MUST BE "LIST".'
          CALL USTOPx(' ')
        END IF
C  Get next file name
        NFILE=1
        GO TO 10
      END IF
C
C8------CHECK FOR "BAS" FILE TYPE.
      IF(FILTYP.EQ.'BAS6') THEN
         INBAS=IU
         FILSTAT='OLD    '
         FILACT=ACTION(1)
C
C9------CHECK FOR "UNFORMATTED" FILE TYPE.
      ELSE IF(FILTYP.EQ.'DATA(BINARY)' .OR.
     1        FILTYP.EQ.'DATAGLO(BINARY)') THEN
         FMTARG=FORM
         ACCARG=ACCESS
C
C10-----CHECK FOR "FORMATTED" FILE TYPE.
      ELSE IF(LINE(ITYP1:ITYP2).EQ.'DATA' .OR.
     1        LINE(ITYP1:ITYP2).EQ.'DATAGLO') THEN
         FMTARG='FORMATTED'
         ACCARG='SEQUENTIAL'
C
C11-----CHECK FOR MAJOR OPTIONS.
      ELSE
        DO 20 I=1,NIUNIT
           IF(LINE(ITYP1:ITYP2).EQ.CUNIT(I)) THEN
              IUNIT(I)=IU
              FILSTAT='OLD    '
              FILACT=ACTION(1)
              GO TO 30
           END IF
20      CONTINUE
        WRITE(IOUT,21) LINE(ITYP1:ITYP2)
21      FORMAT(1X,'ILLEGAL FILE TYPE IN NAME FILE: ',A)
        msg ='Illegal file type in name file (' // LINE(ITYP1:ITYP2) //
     &   ') will be ignored.'
        call store_warning(msg)
30      CONTINUE
        !
        ! Create a package writer if file is input for a package
        select case (FILTYP)
        case ('CHD')
          PkgWriter => model%AddPackage('CHD')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('DRN')
          PkgWriter => model%AddPackage('DRN')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('EVT')
          PkgWriter => model%AddPackage('EVT')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('ETS')
          ! Convert ETS to EVT
          PkgWriter => model%AddPackage('ETS')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('GHB')
          PkgWriter => model%AddPackage('GHB')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('LAK')
          PkgWriter => model%AddPackage('LAK')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('MNW2')
          PkgWriter => model%AddPackage('MAW')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('RCH')
          PkgWriter => model%AddPackage('RCH')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('RIV')
          PkgWriter => model%AddPackage('RIV')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('SFR')
!          ! Convert SFR to SFR6
          PkgWriter => model%AddPackage('SFR')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('UZF')
          PkgWriter => model%AddPackage('UZF')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('WEL')
          PkgWriter => model%AddPackage('WEL')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('HFB6')
          PkgWriter => model%AddPackage('HFB')
          if (associated(PkgWriter)) then
            PkgWriter%IuOrig = iu
          endif
        case ('FHB')
          model%FhbWriter%Active = .true.
          model%FhbWriter%IuOrig = iu
        case ('DE4')
          model%ImsWriter%mf2005_solver = 'DE4'
          model%ImsWriter%IuOrig = iu
        case ('GMG')
          model%ImsWriter%mf2005_solver = 'GMG'
          model%ImsWriter%IuOrig = iu
        case ('PCG')
          model%ImsWriter%mf2005_solver = 'PCG'
          model%ImsWriter%IuOrig = iu
        case ('PCGN')
          model%ImsWriter%mf2005_solver = 'PCGN'
          model%ImsWriter%IuOrig = iu
        case ('SIP')
          model%ImsWriter%mf2005_solver = 'SIP'
          model%ImsWriter%IuOrig = iu
        case ('NWT')
          model%ImsWriter%mf2005_solver = 'NWT'
          model%ImsWriter%IuOrig = iu
        case ('LPF')
          model%NpfWriter%FlowPackage = 'LPF'
        case ('BCF6')
          model%NpfWriter%FlowPackage = 'BCF'
        case ('UPW')
          model%NpfWriter%FlowPackage = 'UPW'
        case ('DIS')
          ! Handled elsewhere
          continue
        case ('PVAL', 'MULT', 'ZONE')
          ! Ignore--parameter values are substituted by RP routines.
          continue
        case ('OC')
          model%OcWriter%Active = .true.
        case ('CHOB')
          if (.not.
     &        associated(model%IbChdWriter%PkgObsWriter)) then
            call createChdObsWriter(newChdObsWriter, model%BaseName, iu)
            model%IbChdWriter%PkgObsWriter => newChdObsWriter
          else
!            call model%IbChdWriter%PkgObsWriter%InitializeObs(
!     &                    model%BaseName)
            model%IbChdWriter%PkgObsWriter%IuObs = iu
          endif
          iuchob = iu
          model%IbChdWriter%PkgObsWriter%Mf6Files => model%Mf6Files
          model%IbChdWriter%PkgObsWriter%Active = .true.
          model%IbChdWriter%ObsActive = .true.
          model%IuChObs = iu
        case ('DROB')
          model%IuDrnObs = iu
        case ('DTOB')
          model%IuDrtObs = iu
        case ('GBOB')
          model%IuGhbObs = iu
        case ('HOB')
          model%HedObsWriter%Active = .true.
          model%IuHeadObs = iu
        case ('RVOB')
          model%IuRivObs = iu
        case default
          msg = 'Conversion of file type ' //
     &          trim(filtyp) // ' is not yet supported.'
          call store_warning(msg)
        end select
      END IF
C
C12-----FOR DATA FILES, CHECK FOR "REPLACE" OR "OLD" OPTION
      IF (FILSTAT.EQ.'UNKNOWN') THEN
        CALL URWORD(LINE,LLOC,IOPT1,IOPT2,1,N,R,IOUT,INUNIT)
        IF (LINE(IOPT1:IOPT2).EQ.'REPLACE' .OR.
     &      LINE(IOPT1:IOPT2).EQ.'OLD')
     &      FILSTAT = LINE(IOPT1:IOPT2)
      ENDIF
C12A----Open file as read-only when 'OLD' is present to allow parallel
C12A----model runs to read data from file simultaneously.
      IF (FILACT.EQ.' ') THEN
        IF (FILSTAT.EQ.'OLD') THEN
          FILACT=ACTION(1)
        ELSE
          FILACT=ACTION(2)
        ENDIF
      ENDIF
C
C13-----WRITE THE FILE NAME AND OPEN IT.
50    FORMAT(1X,/1X,'OPENING ',A,/
     &  1X,'FILE TYPE:',A,'   UNIT ',I4,3X,'STATUS:',A,/
     &  1X,'FORMAT:',A,3X,'ACCESS:',A)
      if (FILSTAT.EQ.'OLD') then
        WRITE(IOUT,50) FNAME(1:IFLEN),
     1     LINE(ITYP1:ITYP2),IU,FILSTAT,FMTARG,ACCARG
        OPEN(UNIT=IU,FILE=FNAME(1:IFLEN),FORM=FMTARG,
     1      ACCESS=ACCARG,STATUS=FILSTAT,ACTION=FILACT,ERR=2000)
        if (FILTYP=='DATA(BINARY)') then
          fcodelocal = FCDATABIN
        else
          fcodelocal = FCDATAIN
        endif
        ! Add input file to list of MF2005 files
        call model%Mf2005Files%AddFile(fname, filtyp, iu, fcodelocal)
      elseif (FILSTAT=='REPLACE') then
        ! Do not open MF2005 output file,
        ! but add it to list of MF6 output files.
        if (FILTYP=='DATA(BINARY)') then
          fcodelocal = FCDATABOUT
        else
          fcodelocal = FCDATAOUT
        endif
        ! Modify file name by prepending base-name
        fname = trim(model%BaseName) // '_' // fname
        call model%Mf6Files%AddFile(FNAME,FILTYP,iu,fcodelocal)
      else
        msg = 'FILSTAT is neither OLD nor REPLACE for file: ' //
     &        trim(fname)
        call store_error(msg)
        call ustop()
      endif
      NFILE=NFILE+1
      GO TO 10
      !
1000  continue
      !
      if (iuchob > 0) then
        ! Go through model packages to see if CHD pkg is active. If so,
        ! set up checks for those CHDs for CHD observation locations.
        npkgs = model%PackageWriters%Count()
        do i=1,npkgs
          PkgWriter => model%GetPackageWriter(i)
          if (PkgWriter%source == 'CHD') then
            select type (PkgWriter)
            type is (ChdPackageWriterType)
              chdpw => PkgWriter
              if (.not.
     &            associated(chdpw%PkgObsWriter)) then
                call createChdObsWriter(newChdObsWriter, model%BaseName,
     &                                  iuchob)
                chdpw%PkgObsWriter => newChdObsWriter
              else
                chdpw%PkgObsWriter%IuObs = iuchob
              endif
              chdpw%PkgObsWriter%Source = 'CHD'
              chdpw%PkgObsWriter%Mf6Files => model%Mf6Files
              chdpw%PkgObsWriter%Active = .true.
              chdpw%ObsActive = .true.
              model%IuChObs = iuchob
            end select
          endif
        enddo
      endif
      !
C
C14-----END OF NAME FILE.  RETURN PROVIDED THAT LISTING FILE AND BAS
C14-----FILES HAVE BEEN OPENED.
      IF(NFILE.EQ.0) THEN
         WRITE(*,*) ' NAME FILE IS EMPTY.'
         CALL USTOPx(' ')
      ELSE IF(INBAS.EQ.0) THEN
         WRITE(IOUT,*) ' BAS PACKAGE FILE HAS NOT BEEN OPENED.'
         CALL USTOPx(' ')
      END IF
      !
      ! If solver has not been defined, use default SMS input
      if (model%ImsWriter%mf2005_solver == '') then
        model%ImsWriter%mf2005_solver = 'DFLT'
        FNAME = trim(model%BaseName) // '.ims6'
        FILTYP = 'IMS6'
        fcodelocal = FCINPUT
        call model%Mf6Files%AddFile(FNAME,FILTYP,iu,fcodelocal)
      endif
      CLOSE (UNIT=INUNIT)
C
      RETURN
C
C15-----ERROR OPENING FILE.
 2000 CONTINUE
      WRITE(*,2010)FNAME(1:IFLEN),IU,FILSTAT,FMTARG,ACCARG,FILACT
      WRITE(IOUT,2010)FNAME(1:IFLEN),IU,FILSTAT,FMTARG,ACCARG,FILACT
 2010 FORMAT(/,1X,'*** ERROR OPENING FILE "',A,'" ON UNIT ',I5,/,
     &7X,'SPECIFIED FILE STATUS: ',A,/
     &7X,'SPECIFIED FILE FORMAT: ',A,/
     &7X,'SPECIFIED FILE ACCESS: ',A,/
     &7X,'SPECIFIED FILE ACTION: ',A,/
     &2X,'-- STOP EXECUTION (SGWF2BAS7OPEN)')
      CALL USTOPx(' ')
C
      END SUBROUTINE SGWF2BAS7OPEN

C###############################################################################

      SUBROUTINE SGWF2BAS7ARMZ(INZONE,INMULT)
C     ******************************************************************
C     ALLOCATE AND READ MULTIPLIER AND ZONE ARRAYS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,     ONLY:NCOL,NROW,IOUT
      USE PARAMMODULE,ONLY:NZONAR,NMLTAR,ZONNAM,MLTNAM,IZON,RMLT
C
      CHARACTER*20 RW
      CHARACTER*1 COP
      CHARACTER*24 ANAME
      CHARACTER*10 CTMP1,CTMP2
      CHARACTER*200 LINE
      double precision :: r
C     ------------------------------------------------------------------
C
C1------Read Number of Zone Arrays if Zone Option is active.
      NZONAR=0
      IF(INZONE.GT.0) THEN
         WRITE(IOUT,1) INZONE
    1    FORMAT(1X,/1X,'ZONE OPTION, INPUT READ FROM UNIT ',I4)
         CALL URDCOM(INZONE,IOUT,LINE)
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NZONAR,R,IOUT,INZONE)
         WRITE(IOUT,2) NZONAR
    2    FORMAT(1X,I5,' ZONE ARRAYS')
         IF(NZONAR.LT.0) NZONAR=0
      END IF
C
C2------Allocate memory for zone arrays.  Allocate one array element if
C2------there are no zone arrays.
      IF(NZONAR.GT.0) THEN
        ALLOCATE (ZONNAM(NZONAR))
        ALLOCATE (IZON(NCOL,NROW,NZONAR))
      ELSE
        ALLOCATE (ZONNAM(1))
        ALLOCATE (IZON(1,1,1))
      ENDIF
C
C3------Read Number of Multiplier Arrays if Multiplier Option is active.
      NMLTAR=0
      IF(INMULT.GT.0) THEN
         WRITE(IOUT,11) INMULT
   11    FORMAT(1X,/1X,'MULTIPLIER OPTION, INPUT READ FROM UNIT ',I4)
         CALL URDCOM(INMULT,IOUT,LINE)
         LLOC=1
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NMLTAR,R,IOUT,INMULT)
         WRITE(IOUT,12) NMLTAR
   12    FORMAT(1X,I3,' MULTIPLIER ARRAYS')
         IF(NMLTAR.LT.0) NMLTAR=0
      END IF
C
C4------Allocate memory for multiplier arrays.  Allocate one array element if
C4------there are no multiplier arrays.
      IF(NMLTAR.GT.0) THEN
        ALLOCATE (MLTNAM(NMLTAR))
        ALLOCATE (RMLT(NCOL,NROW,NMLTAR))
      ELSE
        ALLOCATE (MLTNAM(1))
        ALLOCATE (RMLT(1,1,1))
      ENDIF
C
C5------Initialize names of zones, multipliers, and parameters.
      IF(NZONAR.GT.0) THEN
        DO 10 I=1,NZONAR
        ZONNAM(I)=' '
10      CONTINUE
      END IF
      IF(NMLTAR.GT.0) THEN
        DO 20 I=1,NMLTAR
        MLTNAM(I)=' '
20      CONTINUE
      END IF
C
C6------Define the multiplier arrays.
      IF(NMLTAR.GT.0) THEN
        DO 2000 M=1,NMLTAR
C
C6A-----Read a line describing a multiplier array.
          READ (INMULT,'(A)') LINE
C
C6B-----Get the name of the new array
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,INMULT)
C
C6C-----Add new multiplier name into list.
          MLTNAM(M)=LINE(ISTART:ISTOP)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,INMULT)
          IF(LINE(ISTART:ISTOP).NE.'FUNCTION') THEN
C
C6D-----Define array using array reader.
             ANAME=' MULT. ARRAY: '//MLTNAM(M)
             CALL U2DREL(RMLT(:,:,M),ANAME,NROW,NCOL,0,INMULT,IOUT)
          ELSE
C
C6E-----Define array as aritmetic combination of other multiplier arrays.
C6E-----Start by initializing the array to 0.
             WRITE(IOUT,30) MLTNAM(M)
   30        FORMAT(1X,/1X,'Calculated multiplier array: ',A)
             DO 40 I=1,NROW
             DO 40 J=1,NCOL
             RMLT(J,I,M)=0.
   40        CONTINUE
C
C6E1----Get the names of the multipliers and the operands.
             READ (INMULT,'(A)') LINE
             LLOC=1
             NOP=0
C
C6E2----Get the operator.
   45        IF(NOP.EQ.0) THEN
C
C6E2A---No operator is specified before the first operand -- define it to be " "
                COP=' '
             ELSE
C
C6E2B---Get the operator that precedes each operand after the first operand.
                CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,INMULT)
                IF(LINE(ISTART:ISTOP).EQ.'+' .OR.
     1             LINE(ISTART:ISTOP).EQ.'-' .OR.
     2             LINE(ISTART:ISTOP).EQ.'*' .OR.
     3             LINE(ISTART:ISTOP).EQ.'/') THEN
                   COP=LINE(ISTART:ISTOP)
                ELSE
                   GO TO 1000
                END IF
             END IF
             NOP=NOP+1
C
C6E3----Get the operand.
             CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,INMULT)
             WRITE(IOUT,47 ) COP,LINE(ISTART:ISTOP)
   47        FORMAT(1X,'                        ',A,' ARRAY ',A)
C
C6E4----Lookup the operand in the list of existing multipliers
             CTMP2=LINE(ISTART:ISTOP)
             CALL UPCASE(CTMP2)
             DO 50 MM=1,M-1
               CTMP1=MLTNAM(MM)
               CALL UPCASE(CTMP1)
               IF(CTMP1.EQ.CTMP2) GO TO 60
   50        CONTINUE
             WRITE(IOUT,51) LINE(ISTART:ISTOP)
   51        FORMAT(1X,
     1        'ARRAY OPERAND HAS NOT BEEN PREVIOUSLY DEFINED:',A)
             CALL USTOPx(' ')
C
C6E5----Apply the + operator.
   60        IF(COP.EQ.'+' .OR. COP.EQ.' ') THEN
                DO 100 I = 1, NROW
                DO 100 J = 1, NCOL
                  RMLT(J,I,M) = RMLT(J,I,M)+ RMLT(J,I,MM)
  100           CONTINUE
             ELSE IF(COP.EQ.'-') THEN
                DO 200 I = 1, NROW
                DO 200 J = 1, NCOL
                  RMLT(J,I,M) = RMLT(J,I,M)- RMLT(J,I,MM)
  200           CONTINUE
             ELSE IF(COP.EQ.'*') THEN
                DO 300 I = 1, NROW
                DO 300 J = 1, NCOL
                  RMLT(J,I,M) = RMLT(J,I,M)* RMLT(J,I,MM)
  300           CONTINUE
             ELSE
                DO 400 I = 1, NROW
                DO 400 J = 1, NCOL
                  RMLT(J,I,M) = RMLT(J,I,M)/ RMLT(J,I,MM)
  400           CONTINUE
             END IF
C
C6E6----Get the next operator.
             GO TO 45
C
C6E7-----Done defining the array.  Get the print code and print the array.
1000          IPRN=0
              L=20-ISTOP+ISTART
              IF(L.GT.1)  THEN
                 RW=' '
                 RW(L:20)=LINE(ISTART:ISTOP)
                 READ(RW,'(I20)',ERR=1200) IPRN
              END IF
 1200         IF(IPRN.GE.0) THEN
                 ANAME=' MULT. ARRAY: '//MLTNAM(M)
                 CALL ULAPRWC(RMLT(:,:,M),NCOL,NROW,0,IOUT,IPRN,
     1                 ANAME)
              END IF
          END IF
 2000   CONTINUE
      ENDIF
C
C7------Read the zone array names and arrays
      IF(NZONAR.GT.0) THEN
         DO 3000 NZ=1,NZONAR
         READ(INZONE,'(A)') ZONNAM(NZ)
         CALL U2DINT(IZON(:,:,NZ),'  ZONE ARRAY: '//ZONNAM(NZ),
     1            NROW,NCOL,0,INZONE,IOUT)
 3000    CONTINUE
      END IF
C
C8------Return.
      RETURN
      END SUBROUTINE SGWF2BAS7ARMZ

C###############################################################################

      SUBROUTINE SGWF2BAS7ARPVAL(IUPVAL)
C     ******************************************************************
C     READ PARAMETER INPUT FILE
C     ******************************************************************
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,    ONLY: IOUT,IUNIT
      USE PARAMMODULE, ONLY:MXPAR,IPSUM,PARNAM,B,NPVAL
C
      CHARACTER*10 PNI, PNJ
      CHARACTER*200 LINE
      double precision :: dum
C     ------------------------------------------------------------------
C
C1------CHECK TO SEE IF THE PARAMETER FILE WAS DECLARED IN THE NAME FILE.
      IU=IUNIT(IUPVAL)
      IF(IU.LE.0) THEN
         NPVAL=0
         RETURN
      END IF
C
C2------INITIALIZE VARIABLES
      IERR = 0
      NPE = 0
C
C3------IDENTIFY PARAMETER VALUE OPTION.
      WRITE (IOUT,12) IU
   12 FORMAT (1X,/,1X,
     1  'PARAMETER VALUE INPUT FILE,  INPUT READ FROM UNIT ',I4)
C
C4------READ & PRINT NUMBER OF PARAMETER VALUES.
      CALL URDCOM(IU,IOUT,LINE)
      LLOC = 1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPVAL,DUM,IOUT,IU)
      WRITE (IOUT,14) NPVAL
   14 FORMAT (1X,/,1X,'NUMBER OF PARAMETER VALUES TO BE READ FROM',
     1               ' PARAMETER VALUE FILE:',I5)
      IF (NPVAL.LE.0) THEN
        WRITE (IOUT,16)
   16   FORMAT(1X,'NPVAL IN PARAMETER INPUT FILE MUST BE',
     1         ' > 0 -- STOP EXECUTION')
        CALL USTOPx(' ')
      ENDIF
      IPSUM=NPVAL
C
C5-----DEACTIVATE OPTION IF THERE ARE NO PARAMETERS IN FILE.
      IF(NPVAL.LE.0) THEN
         WRITE(IOUT,*) ' NPVAL in parameter file is 0,',
     1            ' so ignoring the parameter file'
        CLOSE(UNIT=IU)
        IU=0
        RETURN
      END IF
C
C6------STOP IF THERE ARE MORE THAN THE MAXIMUM NUMBER OF PARAMETERS.
      IF(NPVAL.GT.MXPAR) THEN
         WRITE(IOUT,*) ' PARAMETER FILE CONTAINS',NPVAL,
     1     ' VALUES, BUT THE MAXIMUM NUMBER OF PARAMETERS IS',MXPAR
         CALL USTOPx(' ')
      END IF
C
C7------WRITE A HEADING FOR THE LIST OF PARAMETERS.
      WRITE (IOUT,520)
  520 FORMAT (/,' INFORMATION ON PARAMETERS LISTED IN PARAMETER FILE',/,
     &             13X,'  VALUE IN',/,
     &   '    NAME     PARAMETER FILE',/,
     &   ' ----------  --------------')
C
C8-----READ AND WRITE PARAMETER NAMES AND VALUES.
      DO 70 I=1,NPVAL
        READ(IU,*,ERR=80) PARNAM(I),B(I)
        WRITE(IOUT,570) PARNAM(I),B(I)
  570   FORMAT(1X,A10,2X,G12.5)
C
C8A-----CHECK FOR DUPLICATE PARAMETER NAME FOR ALL BUT THE FIRST PARAMETER.
        IF (I.GT.1) THEN
          PNI=PARNAM(I)
          CALL UPCASE(PNI)
          IM1 = I-1
          DO 60 J=1,IM1
            PNJ=PARNAM(J)
            CALL UPCASE(PNJ)
            IF (PNI.EQ.PNJ) THEN
              WRITE(IOUT,500) PARNAM(I)
  500         FORMAT (' PARAMETER "',A10,
     &        '" IS LISTED MORE THAN ONCE IN PARAMETER FILE',/,
     &        ' -- STOP EXECUTION')
                IERR = 1
            ENDIF
   60     CONTINUE
        ENDIF
   70 CONTINUE
C
C9------WRITE A MESSAGE EXPLAINING THAT THE PARAMETER VALUES REPLACE THE
C9------VALUES FROM PACKAGE INPUT FILES..
      WRITE (IOUT,620)
  620 FORMAT(1X,77('-'))
      WRITE (IOUT,630)
  630 FORMAT(' FOR THE PARAMETERS LISTED IN THE TABLE ABOVE,',
     &       ' PARAMETER VALUES IN INDIVIDUAL',/,
     &       ' PACKAGE INPUT FILES ARE REPLACED BY THE VALUES FROM',
     &       ' THE PARAMETER INPUT FILE.')
C
C10-----STOP IF THERE WERE DUPLICATE NAMES.
      IF (IERR.GT.0) THEN
        WRITE(IOUT,680)
  680 FORMAT(/,
     &' ERROR FOUND IN PARAMETER INPUT FILE.  SEARCH ABOVE',/,
     &' FOR "STOP EXECUTION"')
         CALL USTOPx(' ')
      ENDIF
C
C11-----CLOSE FILE AND RETURN.
      CLOSE(UNIT=IU)
      RETURN
C
C
   80 WRITE(IOUT,590)
  590 FORMAT(1X,/,1X,
     1  'ERROR ENCOUNTERED IN READING PARAMETER INPUT FILE',/,
     2       ' -- STOP EXECUTION')
      CALL USTOPx(' ')
C
      END SUBROUTINE SGWF2BAS7ARPVAL

C###############################################################################

      SUBROUTINE SGWF2BAS7STPVAL()
C     ******************************************************************
C     CHECK THAT PARAMETER DEFINITIONS ARE COMPLETE.
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY: IOUT
      USE PARAMMODULE, ONLY:NPVAL,PARTYP,PARNAM
C     ------------------------------------------------------------------
      IF(NPVAL.LE.0) RETURN
      IERR=0
C
Cx------CHECK THAT ALL PARAMETERS IN PARAMETER INPUT FILE HAVE BEEN DEFINED.
      DO 90 IP=1,NPVAL
        IF (PARTYP(IP).EQ.' ') THEN
          IERR = 1
          WRITE(IOUT,110) PARNAM(IP)
  110     FORMAT(1X,/,1X,'PARAMETER "',A10,
     1      '" IN PARAMETER INPUT FILE HAS NOT BEEN DEFINED',/,
     2           ' -- STOP EXECUTION')
        ENDIF
   90 CONTINUE
C
      IF(IERR.NE.0) CALL USTOPx(' ')
C
Cx------RETURN.
      RETURN
      END SUBROUTINE SGWF2BAS7STPVAL

      end module GwfBasSubs
