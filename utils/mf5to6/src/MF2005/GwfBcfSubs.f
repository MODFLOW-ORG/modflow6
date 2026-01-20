      module GwfBcfSubs
        use ConstantsModule, only: DZERO
        use ConstantsPHMFModule, only: HDRYDEFAULT
        use SimPHMFModule, only: ustop, store_error, store_note,
     &                       store_warning
        use utl7module, only: U1DREL, U2DREL
        use ModelModule, only: ModelType
        use GWFBCFMODULE, only: SGWF2BCF7PSV

      contains

      SUBROUTINE GWF2BCF7AR(IN,model)
C     ******************************************************************
C     ALLOCATE ARRAYS AND READ DATA FOR BLOCK-CENTERED FLOW PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,
     1                      CC,CV,IFREFM,botm,lbotm,ibound
      USE GWFBASMODULE,ONLY:HDRY
      USE GWFBCFMODULE,ONLY:IBCFCB,IWDFLG,IWETIT,IHDWET,WETFCT,
     1                      LAYCON,LAYAVG,HY,SC1,SC2,WETDRY,CVWD,TRPY,
     &                      vcont
C
      implicit none
      type(ModelType), pointer, intent(inout) :: model
      integer :: in, inam, iss, kb, kk, kt, l, nbot, ntop
      real :: thk
      CHARACTER*24 ANAME(7)
      CHARACTER*12 AVGNAM(4)
      character(len=300) :: msg
      integer :: i, j, k
      integer :: ianyharmonic, ianylog, ianyarith, ianyarithlog, iavg
      logical :: Note10Written = .false.
      DATA AVGNAM/'HARMONIC    ','ARITHMETIC  ',
     1            'LOGARITHMIC ','*UNCONFINED*'/
C
      DATA ANAME(1) /'    PRIMARY STORAGE COEF'/
      DATA ANAME(2) /'    TRANSMIS. ALONG ROWS'/
      DATA ANAME(3) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(4) /'VERT HYD COND /THICKNESS'/
      DATA ANAME(5) /'  SECONDARY STORAGE COEF'/
      DATA ANAME(6) /'COLUMN TO ROW ANISOTROPY'/
      DATA ANAME(7) /'        WETDRY PARAMETER'/
      ! formats
   10 format('In NPF, storage changes under unconfined conditions ',
     &'account for both specific yield and specific storage.',
     &' As a consequence, numerical results of MODFLOW 6 may',
     &' differ from the results of MODFLOW-2005.')
   20 format('Specific storage (SS) has been assigned',
     &' as zero for',
     &' layer ',i0,' of the MODFLOW 6 model because LAYCON equals',
     &' 1 in the MODFLOW-2005 model.')
C     ------------------------------------------------------------------
C1------ALLOCATE SCALAR VARIABLES IN FORTRAN MODULE.
      ALLOCATE(IBCFCB,IWDFLG,IWETIT,IHDWET)
      ALLOCATE(WETFCT)
      ALLOCATE(LAYCON(NLAY),LAYAVG(NLAY))
      call model%NpfWriter%AllocateArrays()
C
      !
      write(*,*)'Processing BCF package input...'
C2------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'BCF -- BLOCK-CENTERED FLOW PACKAGE, VERSION 7',
     1', 5/2/2005',/,9X,'INPUT READ FROM UNIT',I3)
      !
      if (nlay>1) then
        msg = 'Conversion of BCF is non-unique when NLAY > 1.' //
     &        ' Numeric values are provided in MODFLOW 6 input,' //
     &        ' but some of those values may need to be reviewed.' //
     &        ' See notes above for details.'
        call store_warning(msg)
      endif
      !
      ! -- Assign NPF options corresponding to BCF
      model%NpfWriter%Perched = .true.
      model%NpfWriter%VariableCV = .false.
      model%NpfWriter%Dewatered = .false.
C
C3------READ AND PRINT IBCFCB (FLAG FOR PRINTING
C3------OR UNIT# FOR RECORDING CELL-BY-CELL FLOW TERMS), HDRY
C3------(HEAD AT CELLS THAT CONVERT TO DRY), AND WETTING PARAMETERS.
      IF(IFREFM.EQ.0) THEN
         READ(IN,'(I10,F10.0,I10,F10.0,2I10)')
     1              IBCFCB,HDRY,IWDFLG,WETFCT,IWETIT,IHDWET
      ELSE
         READ(IN,*) IBCFCB,HDRY,IWDFLG,WETFCT,IWETIT,IHDWET
      END IF
      hdry = hdrydefault
      model%NpfWriter%Inpfcb = IBCFCB
      model%NpfWriter%Hdry = HDRY
      if (IWDFLG/=0) model%NpfWriter%Rewet = .true.
      if (WETFCT>0.0d0) model%NpfWriter%Wetfct = WETFCT
      if (IWETIT>1) model%NpfWriter%Iwetit = IWETIT
      if (IHDWET/=0) model%NpfWriter%Ihdwet = IHDWET
C
C3A-----DETERMINE ISS FROM ITRSS
      IF(ITRSS.EQ.0) THEN
         ISS=1
      ELSE
         ISS=0
      END IF
C
C3B-----PRINT VALUES
      IF(ISS.EQ.0) then
        WRITE(IOUT,3)
    3   FORMAT(1X,'TRANSIENT SIMULATION')
        ! For BCF, confined storage is always read as storage coefficient
        model%StoWriter%Isfac = 1 !
      endif
      IF(ISS.NE.0) WRITE(IOUT,4)
    4 FORMAT(1X,'STEADY-STATE SIMULATION')
      IF(IBCFCB.LT.0) WRITE(IOUT,8)
    8 FORMAT(1X,'CONSTANT-HEAD CELL-BY-CELL FLOWS WILL BE PRINTED',
     1     ' WHEN ICBCFL IS NOT 0')
      IF(IBCFCB.GT.0) WRITE(IOUT,9) IBCFCB
    9 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT',I3)
      WRITE(IOUT,11) HDRY
   11 FORMAT(1X,'HEAD AT CELLS THAT CONVERT TO DRY=',G13.5)
      IF(IWDFLG.NE.0) GO TO 35
      WRITE(IOUT,12)
   12 FORMAT(1X,'WETTING CAPABILITY IS NOT ACTIVE')
      GO TO 50
C
   35 WRITE(IOUT,36)
   36 FORMAT(1X,'WETTING CAPABILITY IS ACTIVE')
      IF(IWETIT.LE.0) IWETIT=1
      WRITE(IOUT,37)WETFCT,IWETIT
   37 FORMAT(1X,'WETTING FACTOR=',F10.5,
     1     '     WETTING ITERATION INTERVAL=',I4)
      WRITE(IOUT,38)IHDWET
   38 FORMAT(1X,'FLAG THAT SPECIFIES THE EQUATION TO USE FOR HEAD',
     1    ' AT WETTED CELLS=',I4)
C
C4------READ LAYCON & PRINT TITLE FOR LAYCON TABLE.
   50 IF(IFREFM.EQ.0) THEN
         READ(IN,'(40I2)') (LAYCON(I),I=1,NLAY)
      ELSE
         READ(IN,*) (LAYCON(I),I=1,NLAY)
      END IF
      WRITE(IOUT,52)
   52 FORMAT(1X,5X,'LAYER  LAYER-TYPE CODE     INTERBLOCK T',
     1      /1X,5X,44('-'))
C
C5------LOOP THROUGH LAYERS CALCULATING LAYAVG, PRINTING THE LAYER-TYPE
C5------CODE, AND COUNTING LAYERS THAT NEED TOP & BOT ARRAYS.
      NBOT=0
      NTOP=0
      !
      ianyharmonic = 0
      ianylog = 0
      ianyarith = 0
      ianyarithlog = 0
      !
      DO 100 I=1,NLAY
      IF(LAYCON(I).EQ.30 .OR. LAYCON(I).EQ.32) LAYCON(I)=LAYCON(I)-10
      INAM=LAYCON(I)/10 ! gets left digit of Ltype
      LAYAVG(I)=INAM*10 ! = 0, 10, 20, or 30
      IF(LAYAVG(I).LT.0 .OR. LAYAVG(I).GT.30) THEN
         WRITE(IOUT,53) LAYAVG(I)
   53    FORMAT(1X,'INVALID INTERBLOCK T CODE:',I4)
         CALL USTOP(' ')
      END IF
      LAYCON(I)=LAYCON(I)-LAYAVG(I) ! gets right digit of Ltype
      L=LAYCON(I)
      INAM=INAM+1 ! converts left digit to an index for AVGNAM
      WRITE(IOUT,55) I,L,LAYAVG(I),AVGNAM(INAM)
   55 FORMAT(1X,I9,I13,I11,' -- ',A)
      IF(LAYCON(I).LT.0 .OR. LAYCON(I).GT.3) THEN
         WRITE(IOUT,56) LAYCON(I)
   56    FORMAT(1X,'INVALID LAYER TYPE:',I4)
         CALL USTOP(' ')
      END IF
      !
      ! assign icellavg
      select case (layavg(i))
      case (0)
        ! harmonic mean
        ianyharmonic = 1
      case (10)
        ! arithmetic mean - not supported by MF6
        ianyarith = 1
        msg = 'MODFLOW 6 does not support arithmetic mean option'
     &        // ' for Ltype'
        call store_warning(msg)
      case (20)
        ! logarithmic mean
        ianylog = 1
      case (30)
        ! arithmetic/logarithmic
        ianyarithlog = 1
      case default
        call store_error('Invalid layavg')
        call ustop()
      end select
      !
      ! assign Icelltype
      select case (laycon(i))
      case (0)
        ! confined
        model%NpfWriter%Icelltype(i) = 0   ! saturated thickness held constant
        if (model%StoWriter%Active) then
          model%StoWriter%Iconvert(i) = 0    ! confined storage
        endif
      case (1)
        ! unconfined; T varies, S is constant
        model%NpfWriter%NumConvertible
     &            = model%NpfWriter%NumConvertible + 1
        model%NpfWriter%Icelltype(i) = 1   ! saturated thickness varies
        if (model%StoWriter%Active) then
          model%StoWriter%Iconvert(i) = 1    ! not convertible w/r/t storage; use Sy
        endif
      case (2)
        ! confined/unconfined; T is constant, storage flips b/ Ss and Sy
        ! Icelltype = -1: Saturated thickness varies with computed head unless 
        ! the THICKSTRT option is in effect
        model%NpfWriter%Icelltype(i) = -1  
        model%NpfWriter%ThickStrt = .true. !
        if (model%StoWriter%Active) then
          model%StoWriter%Iconvert(i) = 1    ! storage is convertible
            if (.not. Note10Written) then
              write(msg,10)
              call store_note(msg)
              Note10Written = .true.
            endif
        endif
      case (3)
        ! confined/unconfined: T varies, storage flips b/ Ss and Sy
        model%NpfWriter%NumConvertible
     &            = model%NpfWriter%NumConvertible + 1
        model%NpfWriter%Icelltype(i) = 1   ! convertible
        if (model%StoWriter%Active) then
          model%StoWriter%Iconvert(i) = 1    ! convertible
            if (.not. Note10Written) then
              write(msg,10)
              call store_note(msg)
              Note10Written = .true.
            endif
        endif
      case default
          call store_error('Invalid laycon')
          call ustop()
      end select
C
C5A-----SET GLOBAL HEAD-DEPENDENT THICKNESS FLAGS.
      IF (L.EQ.0) THEN
        LAYHDT(I)=0
        LAYHDS(I)=0
      ELSEIF (L.EQ.1) THEN
        LAYHDT(I)=1
        LAYHDS(I)=0
      ELSEIF (L.EQ.2) THEN
        LAYHDT(I)=0
        LAYHDS(I)=1
      ELSE
        LAYHDT(I)=1
        LAYHDS(I)=1
      ENDIF
C
C5B-----ONLY THE TOP LAYER CAN BE UNCONFINED(LAYCON=1).
      IF(L.NE.1 .OR. I.EQ.1) GO TO 70
      WRITE(IOUT,57)
   57 FORMAT(1X,/1X,'LAYER TYPE 1 IS ONLY ALLOWED IN TOP LAYER')
      CALL USTOP(' ')
C
C5C-----LAYER TYPES 1 AND 3 NEED A BOTTOM. ADD 1 TO KB.
   70 IF(L.EQ.1 .OR. L.EQ.3) NBOT=NBOT+1
C
C5D-----LAYER TYPES 2 AND 3 NEED A TOP. ADD 1 TO KT.
      IF(L.EQ.2 .OR. L.EQ.3) NTOP=NTOP+1
  100 CONTINUE
      !
      ! -- Process flags to determine how
      !    cell_averaging should be assigned.
      iavg = ianyarith + ianyharmonic + ianylog + ianyarithlog
      if (iavg == 1) then
        ! Assignment of averaging method is uniform...which method is used?
        if (ianyarith == 1) then
          msg = 'Arithmethic interblock-transmissivity averaging '
     &         // 'method used in BCF changed to HARMONIC for NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'HARMONIC'
        elseif (ianyharmonic == 1) then
          model%NpfWriter%CellAveraging = 'HARMONIC'
        elseif (ianylog == 1) then
          model%NpfWriter%CellAveraging = 'LOGARITHMIC'
        elseif (ianyarithlog == 1) then
          model%NpfWriter%CellAveraging = 'AMT-LMK'
        endif
      else
        ! Assignment of averaging method is not uniform...choose a
        ! reasonable method based on the methods specified.
        if (ianyharmonic==1) then
          msg = ' Non-uniform interblock-transmissivity averaging '
     &         // 'methods found in BCF input. HARMONIC is used in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'HARMONIC'
        elseif (ianylog == 1) then
          msg = ' Non-uniform interblock-transmissivity averaging '
     &         // 'methods found in BCF input. LOGARITHMIC is used '
     &         // 'in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'LOGARITHMIC'
        elseif (ianyarithlog == 1) then
          msg = ' Non-uniform interblock-transmissivity averaging '
     &         // 'methods found in BCF input. AMT-LMK is used '
     &         // 'in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'AMT-LMK'
        endif
      endif
      !
C
C6------ALLOCATE SPACE FOR ARRAYS.
      IF(ISS.EQ.0) THEN
         ALLOCATE(SC1(NCOL,NROW,NLAY))
      ELSE
         ALLOCATE(SC1(1,1,1))
      END IF
      IF(NTOP.GT.0 .AND. ISS.EQ.0) THEN
         ALLOCATE(SC2(NCOL,NROW,NTOP))
      ELSE
         ALLOCATE(SC2(1,1,1))
      END IF
      ALLOCATE(TRPY(NLAY))
      IF(NBOT.GT.0) THEN
         ALLOCATE(HY(NCOL,NROW,NBOT))
      ELSE
         ALLOCATE(HY(1,1,1))
      END IF
      IF(IWDFLG.NE.0 .AND. NBOT.GT.0) THEN
         ALLOCATE(WETDRY(NCOL,NROW,NBOT))
      ELSE
         ALLOCATE(WETDRY(1,1,1))
      END IF
      IF(IWDFLG.NE.0 .AND. NLAY.GT.1) THEN
         ALLOCATE(CVWD(NCOL,NROW,NLAY-1))
      ELSE
         ALLOCATE(CVWD(1,1,1))
      END IF
      if (nlay > 1) then
        allocate(vcont(ncol,nrow,NLAY-1))
      endif
C
C
C7------READ TRPY
      CALL U1DREL(TRPY,ANAME(6),NLAY,IN,IOUT)
      ! Use TRPY to populate NpfWriter%hani
      do i=1,NLAY
        if (trpy(i) .ne. 1.0d0) then
          model%NpfWriter%UseHani = .true.
          exit
        endif
      enddo
      if (model%NpfWriter%UseHani) then
        allocate(model%NpfWriter%hani(NCOL,NROW,NLAY))
        do k=1,NLAY
          do i=1,NROW
            do j=1,NCOL
              model%NpfWriter%hani(j,i,k) = TRPY(k)
            enddo
          enddo
        enddo
      endif
C
C8------READ ARRAYS FOR EACH LAYER.
      KT=0
      KB=0
      DO 200 K=1,NLAY
      KK=K
C
C8A-----FIND ADDRESS OF EACH LAYER IN THREE DIMENSION ARRAYS.
      IF(LAYCON(K).EQ.1 .OR. LAYCON(K).EQ.3) KB=KB+1
      IF(LAYCON(K).EQ.2 .OR. LAYCON(K).EQ.3) KT=KT+1
C
C8B-----READ PRIMARY STORAGE COEFFICIENT INTO ARRAY SC1 IF TRANSIENT.
      IF(ISS.EQ.0)CALL U2DREL(SC1(:,:,K),ANAME(1),NROW,NCOL,KK,IN,IOUT)
      if (model%StoWriter%Active) then
        if (laycon(k) /= 1) then
          model%StoWriter%Ss(:,:,k) = sc1(:,:,k) ! Primary storage is storage coefficient
        else
          ! laycon=1 means sc1 is always used; store into Sy
          model%StoWriter%Ss(:,:,k) = DZERO      ! Let Ss = zero = Sy for LAYCON = 1
          model%StoWriter%Sy(:,:,k) = sc1(:,:,k) ! Assign Sy = sc1
          write(msg,20)model%NpfWriter%Layptr(k) ! new model layer number
          call store_note(msg)
        endif
      endif
C
C8C-----READ TRANSMISSIVITY INTO ARRAY CC IF LAYER TYPE IS 0 OR 2.
      IF(LAYCON(K).EQ.3 .OR. LAYCON(K).EQ.1) GO TO 105
      CALL U2DREL(CC(:,:,K),ANAME(2),NROW,NCOL,KK,IN,IOUT)
      ! Calculate HK then copy to NpfWriter%HK
1010  format('Thickness <= 0 at (col,row,lay): (',i0,',',i0,',',i0,')')
      do i=1,nrow
        do j=1,ncol
          thk = botm(j,i,LBOTM(k)-1) - botm(j,i,LBOTM(k))
          if (thk>0.0) then
            model%NpfWriter%HK(j,i,k) = cc(j,i,k)/thk
          else
            model%NpfWriter%HK(j,i,k) = 0.0d0
            if (IBOUND(j,i,k)/=0) then
              write(msg,1010)j,i,k
              call store_warning(msg)
            endif
          endif
        enddo
      enddo
      GO TO 110
C
C8D-----READ HYDRAULIC CONDUCTIVITY(HY) IF LAYER TYPE IS 1 OR 3.
  105 CALL U2DREL(HY(:,:,KB),ANAME(3),NROW,NCOL,KK,IN,IOUT)
      ! Copy HY to NpfWriter%HK
      model%NpfWriter%HK(:,:,K) = HY(:,:,KB)
      goto 110
C
C8E-----READ VERTICAL HYCOND/THICK INTO ARRAY CV IF NOT BOTTOM LAYER;
C2E-----MULTIPLIED BY CELL AREA TO CONVERT TO CONDUCTANCE LATER.
! I.e. read vcont
  110 IF(K.EQ.NLAY) GO TO 120
      CALL U2DREL(CV(:,:,K),ANAME(4),NROW,NCOL,KK,IN,IOUT)
      ! Save vcont as read from input
      vcont(:,:,k) = cv(:,:,k)
C
C8F-----READ SECONDARY STORAGE COEFFICIENT INTO ARRAY SC2 IF TRANSIENT
C8F-----AND LAYER TYPE IS 2 OR 3.
  120 IF(LAYCON(K).NE.3 .AND. LAYCON(K).NE.2) GO TO 130
      IF(ISS.EQ.0)CALL U2DREL(SC2(:,:,KT),ANAME(5),NROW,NCOL,KK,IN,IOUT)
      if (model%StoWriter%Active) then
        model%StoWriter%Sy(:,:,k) = sc2(:,:,KT)
      endif
C
C8H-----READ WETDRY CODES IF LAYER TYPE IS 1 OR 3 AND WETTING
C8H-----CAPABILITY HAS BEEN INVOKED (IWDFLG NOT 0).
  130 IF(LAYCON(K).NE.3.AND.LAYCON(K).NE.1)GO TO 200
      IF(IWDFLG.EQ.0)GO TO 200
      CALL U2DREL(WETDRY(:,:,KB),ANAME(7),NROW,NCOL,KK,IN,IOUT)
      ! Copy WETDRY to NpfWriter%WetDry
      !model%NpfWriter%WetDry(:,:,k) = WETDRY(:,:,KB)
  200 CONTINUE
C
C9------PREPARE AND CHECK BCF DATA.
      CALL SGWF2BCF7N(ISS)
C
C10-----SAVE POINTERS FOR GRID AND RETURN.
      CALL SGWF2BCF7PSV(model%IGrid)
      RETURN
      END SUBROUTINE GWF2BCF7AR

C     ******************************************************************

      SUBROUTINE SGWF2BCF7C(K)
C     ******************************************************************
C     COMPUTE BRANCH CONDUCTANCE USING HARMONIC MEAN OF BLOCK
C     CONDUCTANCES -- BLOCK TRANSMISSIVITY IS IN CC UPON ENTRY
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,CR,CC,DELR,DELC
      USE GWFBCFMODULE,ONLY:TRPY
C     ------------------------------------------------------------------
C
      ZERO=0.
      TWO=2.
      YX=TRPY(K)*TWO
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 40 I=1,NROW
      DO 40 J=1,NCOL
      T1=CC(J,I,K)
C
C2------IF T=0 THEN SET CONDUCTANCE EQUAL TO 0. GO ON TO NEXT CELL.
      IF(T1.NE.ZERO) GO TO 10
      CR(J,I,K)=ZERO
      GO TO 40
C
C3------IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
C3------BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
   10 IF(J.EQ.NCOL) GO TO 30
      T2=CC(J+1,I,K)
      CR(J,I,K)=TWO*T2*T1*DELC(I)/(T1*DELR(J+1)+T2*DELR(J))
C
C4------IF THIS IS NOT THE LAST ROW(FRONTMOST) THEN CALCULATE
C4------BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
   30 IF(I.EQ.NROW) GO TO 40
      T2=CC(J,I+1,K)
      CC(J,I,K)=YX*T2*T1*DELR(J)/(T1*DELC(I+1)+T2*DELC(I))
   40 CONTINUE
C
C5------RETURN
      RETURN
      END SUBROUTINE SGWF2BCF7C

C     ******************************************************************

      SUBROUTINE SGWF2BCF7N(ISS)
C     ******************************************************************
C     INITIALIZE AND CHECK BCF DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
!      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,LAYCBD,CC,CV,
!     1                      DELR,DELC,IOUT
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,CC,CV,
     1                      DELR,DELC,IOUT
      USE GWFBASMODULE,ONLY:HNOFLO
      USE GWFBCFMODULE,ONLY:IWDFLG,WETDRY,HY,CVWD,LAYCON,LAYAVG,SC1,SC2
C
      DOUBLE PRECISION HCNV
C     ------------------------------------------------------------------
C
C1------MULTIPLY VERTICAL LEAKANCE BY AREA TO MAKE CONDUCTANCE.
      ZERO=0.
      IF(NLAY.EQ.1) GO TO 20
      K1=NLAY-1
      DO 10 K=1,K1
      DO 10 I=1,NROW
      DO 10 J=1,NCOL
      CV(J,I,K)=CV(J,I,K)*DELR(J)*DELC(I)
   10 CONTINUE
C
C2------IF WETTING CAPABILITY IS ACTIVATED, SAVE CV IN CVWD FOR USE WHEN
C2------WETTING CELLS.
      IF(IWDFLG.EQ.0) GO TO 20
      DO 15 K=1,K1
      DO 15 I=1,NROW
      DO 15 J=1,NCOL
      CVWD(J,I,K)=CV(J,I,K)
   15 CONTINUE
C
C3------IF IBOUND=0, SET CV=0 AND CC=0.
   20 DO 30 K=1,NLAY
      DO 30 I=1,NROW
      DO 30 J=1,NCOL
      IF(IBOUND(J,I,K).NE.0) GO TO 30
      IF(K.NE.NLAY) CV(J,I,K)=ZERO
      IF(K.NE.1) CV(J,I,K-1)=ZERO
      CC(J,I,K)=ZERO
   30 CONTINUE
C
C4------INSURE THAT EACH ACTIVE CELL HAS AT LEAST ONE NON-ZERO
C4------TRANSMISSIVE PROPERTY.
      HCNV=HNOFLO
      KB=0
      DO 60 K=1,NLAY
      IF(LAYCON(K).EQ.1 .OR. LAYCON(K).EQ.3) GO TO 50
C
C4A-----WHEN LAYER TYPE IS 0 OR 2, TRANSMISSIVITY OR CV MUST BE NONZERO.
      DO 45 I=1,NROW
      DO 45 J=1,NCOL
      IF(IBOUND(J,I,K).EQ.0) GO TO 45
      IF(CC(J,I,K).NE.ZERO) GO TO 45
      IF(K.EQ.NLAY) GO TO 41
      IF(CV(J,I,K).NE.ZERO) GO TO 45
   41 IF(K.EQ.1) GO TO 42
      IF(CV(J,I,K-1).NE.ZERO) GO TO 45
   42 IBOUND(J,I,K)=0
      HNEW(J,I,K)=HCNV
      WRITE(IOUT,43) K,I,J
   43 FORMAT(1X,'NODE (LAYER,ROW,COL)',3I4,
     1      ' ELIMINATED BECAUSE ALL CONDUCTANCES TO NODE ARE 0')
   45 CONTINUE
      GO TO 60
C
C4B-----WHEN LAYER TYPE IS 1 OR 3, HY OR CV MUST BE NONZERO.
   50 KB=KB+1
      DO 59 I=1,NROW
      DO 59 J=1,NCOL
C
C4B1----IF WETTING CAPABILITY IS ACTIVE, CHECK CVWD.
      IF(IWDFLG.EQ.0) GO TO 55
      IF(WETDRY(J,I,KB).EQ.ZERO) GO TO 55
      IF(K.EQ.NLAY) GO TO 51
      IF(CVWD(J,I,K).NE.ZERO) GO TO 59
   51 IF(K.EQ.1) GO TO 57
      IF(CVWD(J,I,K-1).NE.ZERO) GO TO 59
      GO TO 57
C
C4B2----WETTING CAPABILITY IS INACTIVE, SO CHECK CV AT ACTIVE CELLS.
   55 IF(IBOUND(J,I,K).EQ.0) GO TO 59
      IF(K.EQ.NLAY) GO TO 56
      IF(CV(J,I,K).NE.ZERO) GO TO 59
   56 IF(K.EQ.1) GO TO 57
      IF(CV(J,I,K-1).NE.ZERO) GO TO 59
C
C4B3----CHECK HYDRAULIC CONDUCTIVITY.
   57 IF(HY(J,I,KB).NE.ZERO) GO TO 59
C
C4B4----HY AND CV ARE ALL 0, SO CONVERT CELL TO NO FLOW.
      IBOUND(J,I,K)=0
      HNEW(J,I,K)=HCNV
      IF(IWDFLG.NE.0) WETDRY(J,I,KB)=ZERO
      WRITE(IOUT,43) K,I,J
   59 CONTINUE
   60 CONTINUE
C
C5------CALCULATE HOR. CONDUCTANCE(CR AND CC) FOR CONSTANT T LAYERS.
      DO 70 K=1,NLAY
      KK=K
      IF(LAYCON(K).EQ.3 .OR. LAYCON(K).EQ.1) GO TO 70
      IF(LAYAVG(K).EQ.0) THEN
         CALL SGWF2BCF7C(KK)
      ELSE IF(LAYAVG(K).EQ.10) THEN
         CALL SGWF2BCF7A(KK)
      ELSE
         CALL SGWF2BCF7L(KK)
      END IF
   70 CONTINUE
C
C6------IF TRANSIENT, LOOP THROUGH LAYERS AND CALCULATE STORAGE
C6------CAPACITY.
      IF(ISS.NE.0) GO TO 100
      KT=0
      DO 90 K=1,NLAY
C
C6A-----MULTIPLY PRIMARY STORAGE COEFFICIENT BY DELR & DELC TO GET
C6A-----PRIMARY STORAGE CAPACITY.
      DO 80 I=1,NROW
      DO 80 J=1,NCOL
      SC1(J,I,K)=SC1(J,I,K)*DELR(J)*DELC(I)
   80 CONTINUE
C
C6B-----IF LAYER IS CONF/UNCONF MULTIPLY SECONDARY STORAGE COEFFICIENT
C6B-----BY DELR AND DELC TO GET SECONDARY STORAGE CAPACITY(SC2).
      IF(LAYCON(K).NE.3 .AND. LAYCON(K).NE.2) GO TO 90
      KT=KT+1
      DO 85 I=1,NROW
      DO 85 J=1,NCOL
      SC2(J,I,KT)=SC2(J,I,KT)*DELR(J)*DELC(I)
   85 CONTINUE
   90 CONTINUE
C
C7------RETURN.
  100 RETURN
      END SUBROUTINE SGWF2BCF7N

C     ******************************************************************

      SUBROUTINE SGWF2BCF7A(K)
C     ******************************************************************
C-------COMPUTE CONDUCTANCE USING ARITHMETIC MEAN TRANSMISSIVITY
C-------ACTIVATED BY LAYAVG=10
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,CR,CC,DELR,DELC
      USE GWFBCFMODULE,ONLY:TRPY
C     ------------------------------------------------------------------
C
      ZERO=0.
      YX=TRPY(K)
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 40 I=1,NROW
      DO 40 J=1,NCOL
      T1=CC(J,I,K)
C
C2------IF T=0 THEN SET CONDUCTANCE EQUAL TO 0. GO ON TO NEXT CELL.
      IF(T1.NE.ZERO) GO TO 10
      CR(J,I,K)=ZERO
      GO TO 40
C
C3------IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
C3------BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
   10 IF(J.EQ.NCOL) GO TO 30
      T2=CC(J+1,I,K)
C3A-----ARITHMETIC MEAN INTERBLOCK TRANSMISSIVITY
      IF(T2.EQ.ZERO) THEN
         CR(J,I,K)=ZERO
      ELSE
         CR(J,I,K)=DELC(I)*(T1+T2)/(DELR(J+1)+DELR(J))
      END IF
C
C4------IF THIS IS NOT THE LAST ROW(FRONTMOST) THEN CALCULATE
C4------BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
   30 IF(I.EQ.NROW) GO TO 40
      T2=CC(J,I+1,K)
      IF(T2.EQ.ZERO) THEN
         CC(J,I,K)=ZERO
      ELSE
         CC(J,I,K)=YX*DELR(J)*(T1+T2)/(DELC(I+1)+DELC(I))
      END IF
   40 CONTINUE
C
C5------RETURN
      RETURN
      END SUBROUTINE SGWF2BCF7A

C     ******************************************************************

      SUBROUTINE SGWF2BCF7L(K)
C     ******************************************************************
C-------COMPUTE CONDUCTANCE USING LOGARITHMIC MEAN TRANSMISSIVITY
C-------ACTIVATED BY LAYAVG=20
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,CR,CC,DELR,DELC
      USE GWFBCFMODULE,ONLY:TRPY
C     ------------------------------------------------------------------
C
      ZERO=0.
      TWO=2.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
      YX=TRPY(K)*TWO
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 40 I=1,NROW
      DO 40 J=1,NCOL
      T1=CC(J,I,K)
C
C2------IF T=0 THEN SET CONDUCTANCE EQUAL TO 0. GO ON TO NEXT CELL.
      IF(T1.NE.ZERO) GO TO 10
      CR(J,I,K)=ZERO
      GO TO 40
C
C3------IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
C3------BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
   10 IF(J.EQ.NCOL) GO TO 30
      T2=CC(J+1,I,K)
      IF(T2.EQ.ZERO) THEN
C3A-----SET TO ZERO AND EXIT IF T2 IS ZERO
         CR(J,I,K)=ZERO
         GO TO 30
      END IF
C3B-----LOGARITHMIC MEAN INTERBLOCK TRANSMISSIVITY
      RATIO=T2/T1
      IF(RATIO.GT.FRAC1.OR.RATIO.LT.FRAC2) THEN
         T=(T2-T1)/LOG(RATIO)
      ELSE
         T=HALF*(T1+T2)
      END IF
      CR(J,I,K)=TWO*DELC(I)*T/(DELR(J+1)+DELR(J))
C
C4------IF THIS IS NOT THE LAST ROW(FRONTMOST) THEN CALCULATE
C4------BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
   30 IF(I.EQ.NROW) GO TO 40
      T2=CC(J,I+1,K)
      IF(T2.EQ.ZERO) THEN
         CC(J,I,K)=ZERO
         GO TO 40
      END IF
      RATIO=T2/T1
      IF(RATIO.GT.FRAC1.OR.RATIO.LT.FRAC2) THEN
         T=(T2-T1)/LOG(RATIO)
      ELSE
         T=HALF*(T1+T2)
      END IF
      CC(J,I,K)=YX*DELR(J)*T/(DELC(I+1)+DELC(I))
   40 CONTINUE
C
C5------RETURN
      RETURN
      END SUBROUTINE SGWF2BCF7L

C     ******************************************************************

      SUBROUTINE SGWF2BCF7U(K)
C     ******************************************************************
C-------COMPUTE CONDUCTANCE USING ARITHMETIC MEAN SATURATED THICKNESS
C-------AND LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
C-------NODE HYDRAULIC CONDUCTIVITY IS IN CC,
C-------NODE SATURATED THICKNESS IS IN BUFF
C-------ACTIVATED BY LAYAVG=30
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,CR,CC,BUFF,DELR,DELC
      USE GWFBCFMODULE,ONLY:TRPY
C     ------------------------------------------------------------------
C
      ZERO=0.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
      YX=TRPY(K)
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 40 I=1,NROW
      DO 40 J=1,NCOL
      T1=CC(J,I,K)
C
C2------IF T=0 THEN SET CONDUCTANCE EQUAL TO 0. GO ON TO NEXT CELL.
      IF(T1.NE.ZERO) GO TO 10
      CR(J,I,K)=ZERO
      GO TO 40
C
C3------IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
C3------BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
   10 IF(J.EQ.NCOL) GO TO 30
      T2=CC(J+1,I,K)
      IF(T2.EQ.ZERO) THEN
C3A-----SET TO ZERO AND EXIT IF T2 IS ZERO
         CR(J,I,K)=ZERO
         GO TO 30
      END IF
C3B-----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
      RATIO=T2/T1
      IF(RATIO.GT.FRAC1.OR.RATIO.LT.FRAC2) THEN
         T=(T2-T1)/LOG(RATIO)
      ELSE
         T=HALF*(T1+T2)
      END IF
C3C-----MULTIPLY LOGARITHMIC K BY ARITHMETIC SAT THICK
      CR(J,I,K)=DELC(I)*T*(BUFF(J,I,K)+BUFF(J+1,I,K))
     *               /(DELR(J+1)+DELR(J))
C
C4------IF THIS IS NOT THE LAST ROW(FRONTMOST) THEN CALCULATE
C4------BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
   30 IF(I.EQ.NROW) GO TO 40
      T2=CC(J,I+1,K)
      IF(T2.EQ.ZERO) THEN
         CC(J,I,K)=ZERO
         GO TO 40
      END IF
      RATIO=T2/T1
      IF(RATIO.GT.FRAC1.OR.RATIO.LT.FRAC2) THEN
         T=(T2-T1)/LOG(RATIO)
      ELSE
         T=HALF*(T1+T2)
      END IF
      CC(J,I,K)=YX*DELR(J)*T*(BUFF(J,I,K)+BUFF(J,I+1,K))
     *            /(DELC(I+1)+DELC(I))
   40 CONTINUE
C
C5------RETURN
      RETURN
      END SUBROUTINE SGWF2BCF7U

C     ******************************************************************


      end module GwfBcfSubs

