      module GwfLpfSubs
        use ConstantsModule, only: DZERO
        use ConstantsPHMFModule, only: HDRYDEFAULT
        use GWFLPFMODULE, only: GWF2LPF7PSV
        use ModelModule, only: ModelType
        use SimPHMFModule, only: ustop, store_error, store_note,
     &                       store_warning
        use utl7module, only: U1DREL, U2DREL, 
     &                        urword, URDCOM

      contains

C*******************************************************************************

      SUBROUTINE GWF2LPF7AR(IN,model)
C     ******************************************************************
C     ALLOCATE AND READ DATA FOR LAYER PROPERTY FLOW PACKAGE
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
!      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,
!     1                      NCNFBD,IBOUND,BUFF,BOTM,NBOTM,DELR,DELC,IOUT
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,
     1                      NCNFBD,IBOUND,BUFF,IOUT
      USE GWFBASMODULE,ONLY:HDRY
      USE GWFLPFMODULE,ONLY:ILPFCB,IWDFLG,IWETIT,IHDWET,
     1                      ISFAC,ICONCV,ITHFLG,NOCVCO,NOVFC,WETFCT,
     2                      LAYTYP,LAYAVG,CHANI,LAYVKA,LAYWET,LAYSTRT,
     3                      LAYFLG,VKA,VKCB,SC1,SC2,HANI,WETDRY,HK
C
      implicit none
      type(ModelType), pointer, intent(inout) :: model
      integer :: in, i, ianame, istart, istop, j, k, khani, kk, lloc,
     &           n, ncnvrt, nhani, nopchk, nphani, nphk, nplpf, npss, 
     &           npsy, npvani, npvk, npvkcb, nwetd
      real :: zero
      CHARACTER*14 LAYPRN(5),AVGNAM(3),TYPNAM(2),VKANAM(2),WETNAM(2),
     1            HANNAM
      DATA AVGNAM/'      HARMONIC','   LOGARITHMIC','     LOG-ARITH'/
      DATA TYPNAM/'      CONFINED','   CONVERTIBLE'/
      DATA VKANAM/'    VERTICAL K','    ANISOTROPY'/
      DATA WETNAM/'  NON-WETTABLE','      WETTABLE'/
      DATA HANNAM/'      VARIABLE'/
      CHARACTER*200 LINE
      CHARACTER*24 ANAME(9),STOTXT
      CHARACTER*4 PTYP
      character(len=1000) :: msg
      integer :: ianyharmonic, ianylog, ianyarithlog, iavg, knew
      double precision :: r
      logical :: NovfcFound, ConstantCvFound, NoCvCorrectionFound
      logical :: Note10Written = .false.
C
      DATA ANAME(1) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(2) /'  HORIZ. ANI. (COL./ROW)'/
      DATA ANAME(3) /'     VERTICAL HYD. COND.'/
      DATA ANAME(4) /' HORIZ. TO VERTICAL ANI.'/
      DATA ANAME(5) /'QUASI3D VERT. HYD. COND.'/
      DATA ANAME(6) /'        SPECIFIC STORAGE'/
      DATA ANAME(7) /'          SPECIFIC YIELD'/
      DATA ANAME(8) /'        WETDRY PARAMETER'/
      DATA ANAME(9) /'     STORAGE COEFFICIENT'/
      ! formats
   10 format('In NPF, storage changes under unconfined conditions ',
     &'account for both specific yield and specific storage.',
     &' As a consequence, numerical results of MODFLOW 6 may',
     &' differ from the results of MODFLOW-2005.')
C     ------------------------------------------------------------------
C1------Allocate scalar data.
      ALLOCATE(ILPFCB,IWDFLG,IWETIT,IHDWET)
      ALLOCATE(ISFAC,ICONCV,ITHFLG,NOCVCO,NOVFC)
      ALLOCATE(WETFCT)
      ZERO=0.
      call model%NpfWriter%AllocateArrays()
C
      !
      write(*,*)'Processing LPF package input...'
C2------IDENTIFY PACKAGE
      WRITE(IOUT,1) IN
    1 FORMAT(1X,/1X,'LPF -- LAYER-PROPERTY FLOW PACKAGE, VERSION 7',
     1', 5/2/2005',/,9X,'INPUT READ FROM UNIT ',I4)
C
C3------READ COMMENTS AND ITEM 1.
      CALL URDCOM(IN,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,ILPFCB,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPLPF,R,IOUT,IN)
C
      hdry = hdrydefault
      model%NpfWriter%Inpfcb = ILPFCB
      model%NpfWriter%Hdry = HDRY
C3A-----WRITE ITEM 1
      IF(ILPFCB.LT.0) WRITE(IOUT,8)
    8 FORMAT(1X,'CONSTANT-HEAD CELL-BY-CELL FLOWS WILL BE PRINTED',
     1  ' WHEN ICBCFL IS NOT 0')
      IF(ILPFCB.GT.0) WRITE(IOUT,9) ILPFCB
    9 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
      WRITE(IOUT,11) HDRY
   11 FORMAT(1X,'HEAD AT CELLS THAT CONVERT TO DRY=',1PG13.5)
      IF(NPLPF.GT.0) THEN
         WRITE(IOUT,15) NPLPF
   15    FORMAT(1X,I5,' Named Parameters     ')
      ELSE
         NPLPF=0
         WRITE(IOUT,'(A)') ' No named parameters'
      END IF
C
C3B-----GET OPTIONS.
      ISFAC=0
      ICONCV=0
      ITHFLG=0
      NOCVCO=0
      NOVFC=0
      NOPCHK=0
      STOTXT=ANAME(6)
      !
      ! -- Record LPF options as specified
      NovfcFound = .false.
      ConstantCvFound = .false.
      NoCvCorrectionFound = .false.
      !
   20 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'STORAGECOEFFICIENT') THEN
         ISFAC=1
         model%StoWriter%Isfac = 1
         STOTXT=ANAME(9)
         WRITE(IOUT,21)
   21    FORMAT(1X,'STORAGECOEFFICIENT OPTION:',/,
     1     1X,'Read storage coefficient rather than specific storage')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'CONSTANTCV') THEN
         ConstantCvFound = .true.
         NovfcFound = .true.
         ICONCV=1
         WRITE(IOUT,23)
   23    FORMAT(1X,'CONSTANTCV OPTION:',/,1X,'Constant vertical',
     1         ' conductance for convertible layers')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'THICKSTRT') THEN
         ITHFLG=1
         WRITE(IOUT,25)
   25    FORMAT(1X,'THICKSTRT OPTION:',/,1X,'Negative LAYTYP indicates',
     1 ' confined layer with thickness computed from STRT-BOT')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOCVCORRECTION') THEN
         NoCvCorrectionFound = .true.
         NOCVCO=1
         WRITE(IOUT,27)
   27    FORMAT(1X,'NOCVCORRECTION OPTION:',/,1X,
     1    'Do not adjust vertical conductance when applying',
     2              ' the vertical flow correction')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOVFC') THEN
         NovfcFound = .true.
         NOVFC=1
         NOCVCO=1
         WRITE(IOUT,29)
   29    FORMAT(1X,'NOVFC OPTION:',/,1X,
     1    'Do not apply the vertical flow correction')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPARCHECK') THEN
         NOPCHK=1
         WRITE(IOUT,30)
   30    FORMAT(1X,'NOPARCHECK  OPTION:',/,1X,
     1    'For data defined by parameters, do not check to see if ',
     2        'parameters define data at all cells')
      END IF
      IF(LLOC.LT.200) GO TO 20
      !
      ! -- Assign NPF options according to user-specified LPF options
      ! -- NPF settings corresponding to LPF defaults:
      model%NpfWriter%Perched = .true.
      model%NpfWriter%VariableCV = .true.
      model%NpfWriter%Dewatered = .true.
      !
      if (NoCvCorrectionFound) then
        model%NpfWriter%VariableCV = .false.
        model%NpfWriter%Dewatered = .false.
      endif
      if (NovfcFound) then
        model%NpfWriter%Perched = .false.
        model%NpfWriter%Dewatered = .false.
        if (ConstantCvFound) then
          model%NpfWriter%VariableCV = .false.
        endif
      endif
C
C4------ALLOCATE AND READ LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET, LAYSTRT.
      ALLOCATE(LAYTYP(NLAY))
      ALLOCATE(LAYAVG(NLAY))
      ALLOCATE(CHANI(NLAY))
      ALLOCATE(LAYVKA(NLAY))
      ALLOCATE(LAYWET(NLAY))
      ALLOCATE(LAYSTRT(NLAY))
      READ(IN,*) (LAYTYP(K),K=1,NLAY)
      READ(IN,*) (LAYAVG(K),K=1,NLAY)
      READ(IN,*) (CHANI(K),K=1,NLAY)
      READ(IN,*) (LAYVKA(K),K=1,NLAY)
      READ(IN,*) (LAYWET(K),K=1,NLAY)
      !
      ! Check CHANI to see if model uses any horizontal anisotropy.
      ! If so, use positive values of CHANI to populate NpfWriter%hani
      do k=1,NLAY
        if (CHANI(k) .ne. 1.0d0) then
          model%NpfWriter%UseHani = .true.
          exit
        endif
      enddo
      if (model%NpfWriter%UseHani) then
        allocate(model%NpfWriter%hani(NCOL,NROW,NLAY))
        do k=1,NLAY
          if (CHANI(k) > 0.0d0) then
            do i=1,NROW
              do j=1,NCOL
                model%NpfWriter%hani(j,i,k) = CHANI(k)
              enddo
            enddo
          endif
        enddo
      endif
      !
      ianyharmonic = 0
      ianylog = 0
      ianyarithlog = 0
      !
      do k=1,nlay
        ! assign icellavg
        select case (layavg(k))
        case (0)
          ! harmonic mean
          ianyharmonic = 1
        case (1)
          ! logarithmic mean
          ianylog = 1
        case (2)
          ! arithmetic/logarithmic
          ianyarithlog = 1
        case default
          call store_error('Invalid layavg')
          call ustop()
        end select
        !
        ! assign Icelltype
        select case (laytyp(k))
        case (0)
          ! confined
          model%NpfWriter%Icelltype(k) = 0  ! confined
          if (model%StoWriter%Active) then
            model%StoWriter%Iconvert(k) = 0   ! confined
          endif
        case (:-1) ! all values < 0
          ! confined/unconfined; T is constant, storage flips b/ Ss and Sy
          model%NpfWriter%Icelltype(k) = -1  ! convertible
          model%NpfWriter%ThickStrt = .true. !
          if (model%StoWriter%Active) then
            model%StoWriter%Iconvert(k) = 1    ! convertible
            if (.not. Note10Written) then
              write(msg,10)
              call store_note(msg)
              Note10Written = .true.
            endif
          endif
        case (1:)  ! all values > 0
          ! confined/unconfined: T varies, storage flips b/ Ss and Sy
          model%NpfWriter%NumConvertible =
     &            model%NpfWriter%NumConvertible + 1
          model%NpfWriter%Icelltype(k) = 1  ! convertible
          if (model%StoWriter%Active) then
            model%StoWriter%Iconvert(k) = 1   ! convertible
            if (.not. Note10Written) then
              write(msg,10)
              call store_note(msg)
              Note10Written = .true.
            endif
          endif
        case default
          call store_error('Invalid laytyp')
          call ustop()
        end select
      enddo
      !
      ! -- Process flags to determine how
      !    cell_averaging should be assigned.
      iavg = ianyharmonic + ianylog + ianyarithlog
      if (iavg == 1) then
        ! Assignment of averaging method is uniform...which method is used?
        if (ianyharmonic == 1) then
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
     &         // 'methods found in LPF input. HARMONIC is used in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'HARMONIC'
        elseif (ianylog == 1) then
          msg = ' Non-uniform interblock-transmissivity averaging '
     &         // 'methods found in LPF input. LOGARITHMIC is used '
     &         // 'in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'LOGARITHMIC'
        elseif (ianyarithlog == 1) then
          msg = ' Non-uniform interblock-transmissivity averaging '
     &         // 'methods found in LPF input. AMT-LMK is used '
     &         // 'in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'AMT-LMK'
        endif
      endif
      !
!
C
C4A-----PRINT A TABLE OF VALUES FOR LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET.
      WRITE(IOUT,47)
   47 FORMAT(1X,/3X,'LAYER FLAGS:',/1X,
     1 'LAYER       LAYTYP        LAYAVG         CHANI ',
     2 '       LAYVKA        LAYWET',/1X,75('-'))
      DO 50 K=1,NLAY
      WRITE(IOUT,48) K,LAYTYP(K),LAYAVG(K),CHANI(K),LAYVKA(K),LAYWET(K)
   48 FORMAT(1X,I4,2I14,1PE14.3,2I14)
C
C4A1----SET GLOBAL HEAD-DEPENDENT TRANSMISSIVITY AND STORAGE FLAGS.
      IF (LAYTYP(K).NE.0) THEN
        LAYHDT(K)=1
        LAYHDS(K)=1
      ELSE
        LAYHDT(K)=0
        LAYHDS(K)=0
      ENDIF
   50 CONTINUE
C
C4A2----SET LAYSTRT AND RESET LAYTYP IF THICKSTRT OPTION IS ACTIVE.
      DO 60 K=1,NLAY
      LAYSTRT(K)=0
      IF(LAYTYP(K).LT.0 .AND. ITHFLG.NE.0) THEN
         LAYSTRT(K)=1
         LAYTYP(K)=0
         LAYHDT(K)=0
         LAYHDS(K)=0
         WRITE(IOUT,57) K
   57    FORMAT(1X,'Layer',I5,
     1  ' is confined because LAYTYP<0 and THICKSTRT option is active')
      END IF
   60 CONTINUE
C
C4B-----BASED ON LAYTYP, LAYAVG, CHANI, LAYWET, COUNT THE NUMBER OF EACH
C4B-----TYPE OF 2-D ARRAY; CHECK VALUES FOR CONSISTENCY; AND SETUP
C4B-----POINTERS IN LAYTYP, CHANI, AND LAYWET FOR CONVENIENT ACCESS
C4B-----TO SC2, HANI, and WETDRY.  PRINT INTERPRETED VALUES OF FLAGS.
      NCNVRT=0
      NHANI=0
      NWETD=0
      WRITE(IOUT,67)
   67 FORMAT(1X,/3X,'INTERPRETATION OF LAYER FLAGS:',/1X,
     1  '                       INTERBLOCK     HORIZONTAL',
     2  '    DATA IN',/1X,
     3  '        LAYER TYPE   TRANSMISSIVITY   ANISOTROPY',
     4  '   ARRAY VKA   WETTABILITY',/1X,
     5  'LAYER      (LAYTYP)      (LAYAVG)       (CHANI)',
     6  '      (LAYVKA)      (LAYWET)',/1X,75('-'))
      DO 100 K=1,NLAY
      IF(LAYTYP(K).NE.0) THEN
         NCNVRT=NCNVRT+1
         LAYTYP(K)=NCNVRT
      END IF
      IF(CHANI(K).LE.ZERO) THEN
         NHANI=NHANI+1
         CHANI(K)=-NHANI
      END IF
      IF(LAYWET(K).NE.0) THEN
         IF(LAYTYP(K).EQ.0) THEN
            WRITE(IOUT,*)
     1          ' LAYWET is not 0 and LAYTYP is 0 for layer:',K
            WRITE(IOUT,*) ' LAYWET must be 0 if LAYTYP is 0'
            CALL USTOP(' ')
         ELSE
            NWETD=NWETD+1
            LAYWET(K)=NWETD
         END IF
      END IF
      IF(LAYAVG(K).LT.0 .OR. LAYAVG(K).GT.2) THEN
         WRITE(IOUT,74) LAYAVG(K)
   74    FORMAT(1X,I8,
     1    ' IS AN INVALID LAYAVG VALUE -- MUST BE 0, 1, or 2')
         CALL USTOP(' ')
      END IF
      LAYPRN(1)=TYPNAM(1)
      IF(LAYTYP(K).NE.0) LAYPRN(1)=TYPNAM(2)
      LAYPRN(2)=AVGNAM(LAYAVG(K)+1)
      IF(CHANI(K).LE.0) THEN
         LAYPRN(3)=HANNAM
      ELSE
         WRITE(LAYPRN(3),'(1PE14.3)') CHANI(K)
      END IF
      LAYPRN(4)=VKANAM(1)
      IF(LAYVKA(K).NE.0) LAYPRN(4)=VKANAM(2)
      LAYPRN(5)=WETNAM(1)
      IF(LAYWET(K).NE.0) LAYPRN(5)=WETNAM(2)
      WRITE(IOUT,78) K,(LAYPRN(I),I=1,5)
   78 FORMAT(1X,I4,5A)
  100 CONTINUE
C
C4C-----PRINT WETTING INFORMATION.
      IF(NWETD.EQ.0) THEN
         WRITE(IOUT,13)
   13    FORMAT(1X,/,1X,'WETTING CAPABILITY IS NOT ACTIVE IN ANY LAYER')
         IWDFLG=0
         model%NpfWriter%Rewet = .false.
      ELSE
         WRITE(IOUT,12) NWETD
   12    FORMAT(1X,/,1X,'WETTING CAPABILITY IS ACTIVE IN',I4,' LAYERS')
         IWDFLG=1
         READ(IN,*) WETFCT,IWETIT,IHDWET
         IF(IWETIT.LE.0) IWETIT=1
         WRITE(IOUT,*) ' WETTING FACTOR=',WETFCT
         WRITE(IOUT,*) ' WETTING ITERATION INTERVAL=',IWETIT
         WRITE(IOUT,*) ' IHDWET=',IHDWET
         model%NpfWriter%Rewet = .true.
         if (WETFCT>0.0d0) model%NpfWriter%Wetfct = WETFCT
         if (IWETIT>1) model%NpfWriter%Iwetit = IWETIT
         if (IHDWET/=0) model%NpfWriter%Ihdwet = IHDWET
      END IF
C
C5------ALLOCATE MEMORY FOR ARRAYS.
      ALLOCATE(LAYFLG(6,NLAY))
      ALLOCATE(HK(NCOL,NROW,NLAY))
      ALLOCATE(VKA(NCOL,NROW,NLAY))
      IF(NCNFBD.GT.0) THEN
         ALLOCATE(VKCB(NCOL,NROW,NCNFBD))
      ELSE
         ALLOCATE(VKCB(1,1,1))
      END IF
      IF(ITRSS.NE.0) THEN
         ALLOCATE(SC1(NCOL,NROW,NLAY))
      ELSE
         ALLOCATE(SC1(1,1,1))
      END IF
      IF(ITRSS.NE.0 .AND. NCNVRT.GT.0) THEN
         ALLOCATE(SC2(NCOL,NROW,NCNVRT))
      ELSE
         ALLOCATE(SC2(1,1,1))
      END IF
      IF(NHANI.GT.0) THEN
         ALLOCATE(HANI(NCOL,NROW,NHANI))
      ELSE
         ALLOCATE(HANI(1,1,1))
      END IF
      IF(NWETD.GT.0) THEN
         ALLOCATE(WETDRY(NCOL,NROW,NWETD))
      ELSE
         ALLOCATE(WETDRY(1,1,1))
      END IF
C
C6------READ PARAMETER DEFINITIONS
      NPHK=0
      NPVKCB=0
      NPVK=0
      NPVANI=0
      NPSS=0
      NPSY=0
      NPHANI=0
      IF(NPLPF.GT.0) THEN
         WRITE(IOUT,115)
  115    FORMAT(/,' PARAMETERS DEFINED IN THE LPF PACKAGE')
         DO 120 K=1,NPLPF
         CALL UPARARRRP(IN,IOUT,N,1,PTYP,1,0,-1)
C   Note that NPHK and the other NP variables in
C   this group are used only as flags, not counts
         IF(PTYP.EQ.'HK') THEN
            NPHK=1
         ELSE IF(PTYP.EQ.'HANI') THEN
C6A-----WHEN A HANI PARAMETER IS USED, THEN ALL HORIZONTAL ANISOTROPY
C6A-----MUST BE DEFINED USING PARAMETERS.  ENSURE THAT ALL CHANI <= 0
            DO 118 I = 1, NLAY
              IF (CHANI(I).GT.0.0) THEN
                WRITE(IOUT,117)
  117           FORMAT(/,
     &' ERROR: WHEN A HANI PARAMETER IS USED, CHANI FOR ALL LAYERS',/,
     &' MUST BE LESS THAN OR EQUAL TO 0.0 -- STOP EXECUTION',
     &' (GWF2LPF7AR)')
                CALL USTOP(' ')
              ENDIF
  118       CONTINUE
            NPHANI=1
         ELSE IF(PTYP.EQ.'VKCB') THEN
            NPVKCB=1
         ELSE IF(PTYP.EQ.'VK') THEN
            NPVK=1
            CALL SGWF2LPF7CK(IOUT,N,'VK  ')
         ELSE IF(PTYP.EQ.'VANI') THEN
            NPVANI=1
            CALL SGWF2LPF7CK(IOUT,N,'VANI')
         ELSE IF(PTYP.EQ.'SS') THEN
            NPSS=1
         ELSE IF(PTYP.EQ.'SY') THEN
            NPSY=1
         ELSE
            WRITE(IOUT,*) ' Invalid parameter type for LPF Package'
            CALL USTOP(' ')
         END IF
  120    CONTINUE
      END IF
C
C7------DEFINE DATA FOR EACH LAYER -- VIA READING OR NAMED PARAMETERS.
      knew = 0
      DO 200 K=1,NLAY
      knew = knew + 1
      KK=K
C
C7A-----DEFINE HORIZONTAL HYDRAULIC CONDUCTIVITY (HK)
      IF(NPHK.EQ.0) THEN
         CALL U2DREL(HK(:,:,KK),ANAME(1),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(1,K)
         WRITE(IOUT,121) ANAME(1),K,LAYFLG(1,K)
  121    FORMAT(1X,/1X,A,' FOR LAYER',I4,
     1   ' WILL BE DEFINED BY PARAMETERS',/1X,'(PRINT FLAG=',I4,')')
         CALL UPARARRSUB1(HK(:,:,KK),NCOL,NROW,KK,'HK',
     1      IOUT,ANAME(1),LAYFLG(1,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,NLAY,
     1      NROW,'HK  ')
      END IF
      !
      ! Copy HK to NpfWriter%HK -- This accounts for parameter values
      ! read from PVAL file, and zone and mult arrays.
      model%NpfWriter%HK(:,:,K) = HK(:,:,KK)
      !
C
C7B-----READ HORIZONTAL ANISOTROPY IF CHANI IS NON-ZERO
      IF(CHANI(K).LE.ZERO) THEN
        KHANI=-CHANI(K)
        IF(NPHANI.EQ.0) THEN
           CALL U2DREL(HANI(:,:,KHANI),ANAME(2),NROW,NCOL,KK,IN,IOUT)
        ELSE
           READ(IN,*) LAYFLG(6,K)
           WRITE(IOUT,121) ANAME(2),K,LAYFLG(6,K)
           CALL UPARARRSUB1(HANI(:,:,KHANI),NCOL,NROW,KK,'HANI',
     1      IOUT,ANAME(2),LAYFLG(6,KK))
           IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,
     1      NLAY,NROW,'HANI')
        END IF
        ! Populate model%NpfWriter%hani
        if (model%NpfWriter%UseHani) then
          do i=1,NROW
            do j=1,NCOL
              model%NpfWriter%hani(j,i,k) = (hani(j,i,khani))
            enddo
          enddo
        endif
      END IF
C
C7C-----DEFINE VERTICAL HYDRAULIC CONDUCTIVITY OR HORIZONTAL TO VERTICAL
C7C-----ANISOTROPY (VKA).
      IANAME=3
      PTYP='VK'
      IF(LAYVKA(K).NE.0) THEN
         IANAME=4
         PTYP='VANI'
      END IF
      IF(NPVK.EQ.0 .AND. NPVANI.EQ.0) THEN
         CALL U2DREL(VKA(:,:,KK),ANAME(IANAME),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(2,K)
         WRITE(IOUT,121) ANAME(IANAME),K,LAYFLG(2,K)
         CALL UPARARRSUB1(VKA(:,:,KK),NCOL,NROW,KK,PTYP,IOUT,
     1                       ANAME(IANAME),LAYFLG(2,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,NLAY,
     1                       NROW,PTYP)
      END IF
      !
      ! Define VK for mf6
      if (layvka(k) == 0) then
        ! VKA contains vertical hydraulic conductivity
        model%NpfWriter%vk(:,:,knew) = vka(:,:,k)
      else
        ! VKA contains vertical anisotropy. Compute VK from HK and anisotropy.
        ! Anisotropy is ratio HK/VK.
        do i=1,NROW
          do j=1,NCOL
            if (vka(j,i,k) > 0.) then
              r = hk(j,i,k) / vka(j,i,k)
            else
              r = 0.
            end if
            model%NpfWriter%vk(j,i,knew) = r
          enddo
        enddo
      endif
C
C7D-----DEFINE SPECIFIC STORAGE OR STORAGE COEFFICIENT IN ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0) THEN
         IF(NPSS.EQ.0) THEN
            CALL U2DREL(SC1(:,:,KK),STOTXT,NROW,NCOL,KK,IN,IOUT)
         ELSE
            READ(IN,*) LAYFLG(3,K)
            WRITE(IOUT,121) STOTXT,K,LAYFLG(3,K)
            CALL UPARARRSUB1(SC1(:,:,KK),NCOL,NROW,KK,'SS',
     1           IOUT,STOTXT,LAYFLG(3,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,
     1           NLAY,NROW,'SS  ')
         END IF
         ! For current layer, define Ss array of NpfWriter
         ! here (before Ss is multiplied by cell thickness
         ! and SC1 is multiplied by cell area).
         model%StoWriter%Ss(:,:,k) = sc1(:,:,k)
         IF(ISFAC.EQ.0) THEN ! ISFAC=0 indicates Ss is read
            CALL SGWF2LPF7SC(SC1(:,:,KK),KK,1)
         ELSE                ! ISFAC=1 indicates S is read
            CALL SGWF2LPF7SC(SC1(:,:,KK),KK,0)
         END IF
      END IF
C
C7E-----DEFINE SPECIFIC YIELD IN ARRAY SC2 IF TRANSIENT AND LAYER IS
C7E-----IS CONVERTIBLE.
      IF(LAYTYP(K).NE.0) THEN
         IF(ITRSS.NE.0) THEN
            IF(NPSY.EQ.0) THEN
               CALL U2DREL(SC2(:,:,LAYTYP(K)),ANAME(7),NROW,NCOL,KK,IN,
     1                 IOUT)
            ELSE
               READ(IN,*) LAYFLG(4,K)
               WRITE(IOUT,121) ANAME(7),K,LAYFLG(4,K)
               CALL UPARARRSUB1(SC2(:,:,LAYTYP(K)),NCOL,
     1         NROW,KK,'SY',IOUT,ANAME(7),LAYFLG(4,KK))
               IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,
     1           NCOL,NLAY,NROW,'SY  ')
            END IF
            ! For current layer, define Sy array of NpfWriter
            ! here (before Sy is multiplied by cell area).
            model%StoWriter%Sy(:,:,k) = sc2(:,:,laytyp(k))
            CALL SGWF2LPF7SC(SC2(:,:,LAYTYP(K)),KK,0)
         END IF
      END IF
C
C7F-----READ CONFINING BED VERTICAL HYDRAULIC CONDUCTIVITY (VKCB)
      IF(LAYCBD(K).NE.0) THEN
         IF(NPVKCB.EQ.0) THEN
            CALL U2DREL(VKCB(:,:,LAYCBD(K)),ANAME(5),NROW,NCOL,KK,IN,
     1             IOUT)
         ELSE
            READ(IN,*) LAYFLG(5,K)
            WRITE(IOUT,121) ANAME(5),K,LAYFLG(5,K)
            CALL UPARARRSUB1(VKCB(:,:,LAYCBD(K)),NCOL,NROW,KK,
     1         'VKCB',IOUT,ANAME(5),LAYFLG(5,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,
     1         NLAY,NROW,'VKCB')
         END IF
         ! Layer K is underlain by a quasi-3D unit.
         ! Increment knew and populate model%NpfWriter%vk from vkcb
         knew = knew + 1
         do i=1,NROW
           do j=1,NCOL
             model%NpfWriter%vk(j,i,knew) = vkcb(j,i,laycbd(k))
           enddo
         enddo
      END IF
C
C7G-----READ WETDRY CODES IF WETTING CAPABILITY HAS BEEN INVOKED
C7G-----(LAYWET NOT 0).
      IF(LAYWET(K).NE.0) THEN
         CALL U2DREL(WETDRY(:,:,LAYWET(K)),ANAME(8),NROW,NCOL,KK,IN,
     1            IOUT)
      END IF
  200 CONTINUE
C
C8------PREPARE AND CHECK LPF DATA.
      CALL SGWF2LPF7N()
C
C9------RETURN
      CALL GWF2LPF7PSV(model%IGrid)
      RETURN
      END SUBROUTINE GWF2LPF7AR

!***********************************************************************

      SUBROUTINE SGWF2LPF7N()
C     ******************************************************************
C     INITIALIZE AND CHECK LPF DATA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
!      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,LAYCBD,CV,
!     1                      BOTM,NBOTM,DELR,DELC,IOUT
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,LAYCBD,IOUT
      USE GWFBASMODULE,ONLY:HNOFLO
      USE GWFLPFMODULE,ONLY:LAYWET,WETDRY,HK,VKCB,LAYTYP,VKA
C     ------------------------------------------------------------------
C
C1------DEFINE CONSTANTS.
      ZERO=0.
      HCNV=HNOFLO
C
C2------INSURE THAT EACH ACTIVE CELL HAS AT LEAST ONE NON-ZERO
C2------TRANSMISSIVE PARAMETER.
      DO 60 K=1,NLAY
      IF(LAYWET(K).NE.0) THEN
C
C3------WETTING IS ACTIVE.
         DO 40 I=1,NROW
         DO 40 J=1,NCOL
         IF(IBOUND(J,I,K).EQ.0 .AND. WETDRY(J,I,LAYWET(K)).EQ.ZERO)
     1                GO TO 40
C
C3A-----CHECK HORIZONTAL HYDRAULIC CONDUCTIVITY (HK).
         IF(HK(J,I,K).NE.ZERO) GO TO 40
C
C3B-----CHECK VERTICAL HYDRAULIC CONDUCTIVITY AND CONFINING BED
C3B-----VERTICAL HYDRAULIC CONDUCTIVITY.
         IF(NLAY.GT.1) THEN

            IF(VKA(J,I,K).NE.ZERO) THEN
               IF(K.NE.NLAY) THEN
                  IF (VKA(J,I,K+1).NE.ZERO) THEN
                     IF(LAYCBD(K).NE.0) THEN
                        IF(VKCB(J,I,LAYCBD(K)).NE.ZERO) GO TO 40
                     ELSE
                        GO TO 40
                     END IF
                  END IF
               END IF
               IF(K.NE.1) THEN
                  IF (VKA(J,I,K-1).NE.ZERO) THEN
                     IF (LAYCBD(K-1).NE.0) THEN
                        IF(VKCB(J,I,LAYCBD(K-1)).NE.ZERO) GO TO 40
                     ELSE
                        GO TO 40
                     END IF
                  ENDIF
               END IF
            END IF
         END IF
C
C3C-----ALL TRANSMISSIVE TERMS ARE ALL 0, SO CONVERT CELL TO NO FLOW.
         IBOUND(J,I,K)=0
         HNEW(J,I,K)=HCNV
         WETDRY(J,I,LAYWET(K))=ZERO
         WRITE(IOUT,43) K,I,J
   40    CONTINUE
C
      ELSE
C
C4------WETTING IS INACTIVE
         DO 50 I=1,NROW
         DO 50 J=1,NCOL
         IF(IBOUND(J,I,K).EQ.0) GO TO 50
C
C4A-----CHECK HORIZONTAL HYDRAULIC CONDUCTIVITY (HK).
         IF(HK(J,I,K).NE.ZERO) GO TO 50
C
C4B-----CHECK VERTICAL HYDRAULIC CONDUCTIVITY AND CONFINING BED
C4B-----VERTICAL HYDRAULIC CONDUCTIVITY.
         IF(NLAY.GT.1) THEN
            IF(VKA(J,I,K).NE.ZERO) THEN
               IF(K.NE.NLAY) THEN
                  IF (VKA(J,I,K+1).NE.ZERO) THEN
                     IF(LAYCBD(K).NE.0) THEN
                        IF(VKCB(J,I,LAYCBD(K)).NE.ZERO) GO TO 50
                     ELSE
                        GO TO 50
                     END IF
                  END IF
               END IF
               IF(K.NE.1) THEN
                  IF (VKA(J,I,K-1).NE.ZERO) THEN
                     IF (LAYCBD(K-1).NE.0) THEN
                        IF(VKCB(J,I,LAYCBD(K-1)).NE.ZERO) GO TO 50
                     ELSE
                        GO TO 50
                     END IF
                  ENDIF
               END IF
            END IF
         END IF
C
C4C-----ALL TRANSMISSIVE TERMS ARE 0, SO CONVERT CELL TO NO FLOW.
         IBOUND(J,I,K)=0
         HNEW(J,I,K)=HCNV
         WRITE(IOUT,43) K,I,J
   43    FORMAT(1X,'NODE (LAYER,ROW,COL) ',I3,2(1X,I5),
     1 ' ELIMINATED BECAUSE ALL HYDRAULIC',/,
     2 ' CONDUCTIVITIES TO NODE ARE 0')
   50    CONTINUE
      END IF
   60 CONTINUE
C
C5------CALCULATE HOR. CONDUCTANCE(CR AND CC) FOR CONSTANT T LAYERS.
      DO 70 K=1,NLAY
      KK=K
      IF(LAYTYP(K).EQ.0) CALL SGWF2LPF7HCOND(KK,0,0,0)
   70 CONTINUE
C
C6------CALCULATE VERTICAL CONDUCTANCE BETWEEN CONFINED LAYERS.
      IF(NLAY.GT.1) THEN
         DO 10 K=1,NLAY-1
         KK=K
         IF(LAYTYP(K).EQ.0 .AND. LAYTYP(K+1).EQ.0)
     1      CALL SGWF2LPF7VCOND(KK)
   10    CONTINUE
      END IF
C
C7------RETURN.
      RETURN
      END SUBROUTINE SGWF2LPF7N

!***********************************************************************

      SUBROUTINE SGWF2LPF7SC(SC,K,ISPST)
C     ******************************************************************
C     COMPUTE STORAGE CAPACITY
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
!      USE GLOBAL,        ONLY:NCOL,NROW,DELR,DELC,BOTM,LBOTM,LAYCBD
      USE GLOBAL,        ONLY:NCOL,NROW,DELR,DELC,BOTM,LBOTM
C
      double precision ::  SC(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------MULTIPLY SPECIFIC STORAGE BY THICKNESS, DELR, AND DELC TO GET
C1------CONFINED STORAGE CAPACITY.
      IF(ISPST.NE.0) THEN
         DO 80 I=1,NROW
         DO 80 J=1,NCOL
         THICK=BOTM(J,I,LBOTM(K)-1)-BOTM(J,I,LBOTM(K))
         SC(J,I)=SC(J,I)*THICK*DELR(J)*DELC(I)
   80    CONTINUE
      ELSE
C
C2------MULTIPLY SPECIFIC YIELD BY DELR AND DELC TO GET UNCONFINED
C2------STORAGEE CAPACITY(SC2).
         DO 85 I=1,NROW
         DO 85 J=1,NCOL
         SC(J,I)=SC(J,I)*DELR(J)*DELC(I)
   85    CONTINUE
      END IF
C
      RETURN
      END SUBROUTINE SGWF2LPF7SC

!***********************************************************************

      SUBROUTINE SGWF2LPF7HCOND(K,KITER,KSTP,KPER)
C     ******************************************************************
C     COMPUTE HORIZONTAL BRANCH CONDUCTANCE FOR ONE LAYER.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
!      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,IBOUND,HNEW,BOTM,NBOTM,
!     1                      LBOTM,CC,STRT
      USE GLOBAL,      ONLY:IOUT,NCOL,NROW,IBOUND,HNEW,BOTM,
     1                      LBOTM,CC,STRT
      USE GWFBASMODULE,ONLY:HDRY
      USE GWFLPFMODULE,ONLY:LAYWET,IWETIT,LAYTYP,LAYAVG,LAYSTRT
C
      CHARACTER*3 ACNVRT
      DIMENSION ICNVRT(5),JCNVRT(5),ACNVRT(5)
C
C     ------------------------------------------------------------------
C1------INITIALIZE DATA.
      ZERO=0.
      NCNVRT=0
      IHDCNV=0
C
C2------IF LAYER IS WETTABLE CONVERT DRY CELLS TO WET WHEN APPROPRIATE.
      ITFLG=1
      IF(LAYWET(K).NE.0) ITFLG=MOD(KITER,IWETIT)
      IF(ITFLG.EQ.0) CALL SGWF2LPF7WET(K,KITER,KSTP,KPER,
     2      IHDCNV,NCNVRT,ICNVRT,JCNVRT,ACNVRT)
C
C3------LOOP THROUGH EACH CELL, AND CALCULATE SATURATED THICKNESS.
      DO 200 I=1,NROW
      DO 200 J=1,NCOL
C
C3A-----SET STAURATED THICKNESS=0. FOR DRY CELLS.
      IF(IBOUND(J,I,K).EQ.0) THEN
         CC(J,I,K)=ZERO
      ELSE
C
C3B-----CALCULATE SATURATED THICKNESS FOR A WET CELL.
         BBOT=BOTM(J,I,LBOTM(K))
         IF(LAYSTRT(K).NE.0) THEN
            TTOP=STRT(J,I,K)
            IF(BBOT.GT.TTOP) THEN
               WRITE(IOUT,33) K,I,J
   33       FORMAT(1X,/1X,'Negative cell thickness at (layer,row,col)',
     1         I4,',',I5,',',I5)
               WRITE(IOUT,34) TTOP,BBOT
   34          FORMAT(1X,'Initial head, bottom elevation:',1P,2G13.5)
               CALL USTOP(' ')
            END IF
         ELSE
            TTOP=BOTM(J,I,LBOTM(K)-1)
            IF(BBOT.GT.TTOP) THEN
               WRITE(IOUT,35) K,I,J
   35       FORMAT(1X,/1X,'Negative cell thickness at (layer,row,col)',
     1         I4,',',I5,',',I5)
               WRITE(IOUT,36) TTOP,BBOT
   36          FORMAT(1X,'Top elevation, bottom elevation:',1P,2G13.5)
               CALL USTOP(' ')
            END IF
         END IF
         IF(LAYTYP(K).NE.0) THEN
            HHD=HNEW(J,I,K)
            IF(HHD.LT.TTOP) TTOP=HHD
         END IF
         THCK=TTOP-BBOT
         CC(J,I,K)=THCK
C
C
C3C-----WHEN SATURATED THICKNESS <= 0, PRINT A MESSAGE AND SET
C3C-----HNEW=HDRY, SATURATED THICKNESS=0.0, AND IBOUND=0.
         IF(THCK.LE.ZERO) THEN
            CALL SGWF2LPF7WDMSG(1,NCNVRT,ICNVRT,JCNVRT,ACNVRT,IHDCNV,
     1             IOUT,KITER,J,I,K,KSTP,KPER,NCOL,NROW)
            HNEW(J,I,K)=HDRY
            CC(J,I,K)=ZERO
            IF(IBOUND(J,I,K).LT.0) THEN
               WRITE(IOUT,151)
  151          FORMAT(1X,/1X,'CONSTANT-HEAD CELL WENT DRY',
     1          ' -- SIMULATION ABORTED')
                    WRITE(IOUT,*) TTOP, BBOT, THCK
               WRITE(IOUT,152) K,I,J,KITER,KSTP,KPER
  152          FORMAT(1X,'LAYER=',I3,' ROW=',I5,' COLUMN=',I5,
     1    ' ITERATION=',I3,' TIME STEP=',I3,' STRESS PERIOD=',I4)
               CALL USTOP(' ')
            END IF
            IBOUND(J,I,K)=0
         END IF
      END IF
  200 CONTINUE
C
C4------PRINT ANY REMAINING CELL CONVERSIONS NOT YET PRINTED.
      CALL SGWF2LPF7WDMSG(0,NCNVRT,ICNVRT,JCNVRT,ACNVRT,IHDCNV,
     1             IOUT,KITER,J,I,K,KSTP,KPER,NCOL,NROW)
C
C5------CHANGE IBOUND VALUE FOR CELLS THAT CONVERTED TO WET THIS
C5------ITERATION FROM 30000 to 1.
      IF(LAYWET(K).NE.0) THEN
         DO 205 I=1,NROW
         DO 205 J=1,NCOL
         IF(IBOUND(J,I,K).EQ.30000) IBOUND(J,I,K)=1
  205    CONTINUE
      END IF
C
C6------COMPUTE HORIZONTAL BRANCH CONDUCTANCES FROM CELL HYDRAULIC
C6------CONDUCTIVITY, SATURATED THICKNESS, AND GRID DIMENSIONS.
      IF(LAYAVG(K).EQ.0) THEN
         CALL SGWF2LPF7HHARM(K)
      ELSE IF(LAYAVG(K).EQ.1) THEN
         CALL SGWF2LPF7HLOG(K)
      ELSE
         CALL SGWF2LPF7HUNCNF(K)
      END IF
C
C7------RETURN.
      RETURN
      END SUBROUTINE SGWF2LPF7HCOND

!***********************************************************************

      SUBROUTINE SGWF2LPF7WET(K,KITER,KSTP,KPER,IHDCNV,NCNVRT,
     1                        ICNVRT,JCNVRT,ACNVRT)
C
C     ******************************************************************
C     CONVERT DRY CELLS TO WET.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:IOUT,NCOL,NROW,NLAY,HNEW,IBOUND,BOTM,LBOTM
!      USE GWFLPFMODULE,  ONLY:LAYTYP,CHANI,LAYVKA,LAYWET,WETDRY,
!     1                        WETFCT,IHDWET
      USE GWFLPFMODULE,  ONLY:LAYWET,WETDRY,WETFCT,IHDWET
C
      CHARACTER*3 ACNVRT
      DIMENSION ICNVRT(5),JCNVRT(5),ACNVRT(5)
C     ------------------------------------------------------------------
C
C1------LOOP THROUGH ALL CELLS.
      ZERO=0.0
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
C
C2------IF CELL IS DRY AND IF IT IS WETTABLE, CONTINUE CHECKING TO SEE
C2------IF IT SHOULD BECOME WET.
      IF(IBOUND(J,I,K).EQ.0 .AND. WETDRY(J,I,LAYWET(K)).NE.ZERO) THEN
C
C3------CALCULATE WETTING ELEVATION.
         WD=WETDRY(J,I,LAYWET(K))
         IF(WD.LT.ZERO) WD=-WD
         TURNON=BOTM(J,I,LBOTM(K))+WD
C
C4------CHECK HEAD IN CELL BELOW TO SEE IF WETTING ELEVATION HAS BEEN
C4------REACHED.
         IF(K.NE.NLAY) THEN
            HTMP=HNEW(J,I,K+1)
            IF(IBOUND(J,I,K+1).GT.0 .AND. HTMP.GE.TURNON) GO TO 50
         END IF
C
C5------CHECK HEAD IN ADJACENT HORIZONTAL CELLS TO SEE IF WETTING
C5------ELEVATION HAS BEEN REACHED.
         IF(WETDRY(J,I,LAYWET(K)).GT.ZERO) THEN
            IF(J.NE.1) THEN
               HTMP=HNEW(J-1,I,K)
               IF(IBOUND(J-1,I,K).GT.0 .AND. IBOUND(J-1,I,K).NE.30000.
     1                       AND. HTMP.GE.TURNON) GO TO 50
            END IF
            IF(J.NE.NCOL) THEN
               HTMP=HNEW(J+1,I,K)
               IF(IBOUND(J+1,I,K).GT.0 .AND. HTMP.GE.TURNON) GO TO 50
            END IF
            IF(I.NE.1) THEN
               HTMP=HNEW(J,I-1,K)
               IF(IBOUND(J,I-1,K).GT.0 .AND. IBOUND(J,I-1,K).NE.30000.
     1                       AND. HTMP.GE.TURNON) GO TO 50
            END IF
            IF(I.NE.NROW) THEN
               HTMP=HNEW(J,I+1,K)
               IF(IBOUND(J,I+1,K).GT.0 .AND. HTMP.GE.TURNON) GO TO 50
            END IF
         END IF
C
C6------WETTING ELEVATION HAS NOT BEEN REACHED, SO CELL REMAINS DRY.
         GO TO 100
C
C7------CELL BECOMES WET.  PRINT MESSAGE, SET INITIAL HEAD, AND SET
C7------IBOUND.
   50    CALL SGWF2LPF7WDMSG(2,NCNVRT,ICNVRT,JCNVRT,ACNVRT,IHDCNV,
     1             IOUT,KITER,J,I,K,KSTP,KPER,NCOL,NROW)
C
C7A-----USE EQUATION 3A IF IHDWET=0; USE EQUATION 3B IF IHDWET IS NOT 0.
         IF(IHDWET.EQ.0) THEN
            HNEW(J,I,K)=BOTM(J,I,LBOTM(K))+
     1                        WETFCT*(HTMP-BOTM(J,I,LBOTM(K)))
         ELSE
            HNEW(J,I,K)=BOTM(J,I,LBOTM(K))+WETFCT*WD
         END IF
         IBOUND(J,I,K)=30000
      END IF
C
C8------END OF LOOP FOR ALL CELLS IN LAYER.
  100 CONTINUE
C
C9------RETURN.
      RETURN
      END SUBROUTINE SGWF2LPF7WET

!***********************************************************************

      SUBROUTINE SGWF2LPF7WDMSG(ICODE,NCNVRT,ICNVRT,JCNVRT,ACNVRT,
     1             IHDCNV,IOUT,KITER,J,I,K,KSTP,KPER,NCOL,NROW)
C     ******************************************************************
C     PRINT MESSAGE WHEN CELLS CONVERT BETWEEN WET AND DRY.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*3 ACNVRT
      DIMENSION ICNVRT(5),JCNVRT(5),ACNVRT(5)
C     ------------------------------------------------------------------
C
C1------KEEP TRACK OF CELL CONVERSIONS.
      IF(ICODE.GT.0) THEN
         NCNVRT=NCNVRT+1
         ICNVRT(NCNVRT)=I
         JCNVRT(NCNVRT)=J
         IF(ICODE.EQ.1) THEN
            ACNVRT(NCNVRT)='DRY'
         ELSE
            ACNVRT(NCNVRT)='WET'
         END IF
      END IF
C
C2------PRINT A LINE OF DATA IF 5 CONVERSIONS HAVE OCCURRED OR IF ICODE
C2------INDICATES THAT A PARTIAL LINE SHOULD BE PRINTED.
      IF(NCNVRT.EQ.5 .OR. (ICODE.EQ.0 .AND. NCNVRT.GT.0)) THEN
         IF(IHDCNV.EQ.0) WRITE(IOUT,17) KITER,K,KSTP,KPER
   17    FORMAT(1X,/1X,'CELL CONVERSIONS FOR ITER.=',I3,'  LAYER=',
     1       I3,'  STEP=',I3,'  PERIOD=',I4,'   (ROW,COL)')
         IHDCNV=1
         IF (NROW.LE.999 .AND. NCOL.LE.999) THEN
           WRITE(IOUT,18) (ACNVRT(L),ICNVRT(L),JCNVRT(L),L=1,NCNVRT)
   18      FORMAT(1X,3X,5(A,'(',I3,',',I3,')   '))
         ELSE
           WRITE(IOUT,19) (ACNVRT(L),ICNVRT(L),JCNVRT(L),L=1,NCNVRT)
   19      FORMAT(1X,2X,5(A,'(',I5,',',I5,')'))
         ENDIF
         NCNVRT=0
      END IF
C
C3------RETURN.
      RETURN
      END SUBROUTINE SGWF2LPF7WDMSG

!***********************************************************************

      SUBROUTINE SGWF2LPF7HHARM(K)
C     ******************************************************************
C     COMPUTE HORIZONTAL BRANCH CONDUCTANCE USING HARMONIC MEAN OF BLOCK
C     CONDUCTANCES (DISTANCE WEIGHTED HARMONIC MEAN OF TRANSMISSIVITY).
C     CELL THICKNESS IS IN CC UPON ENTRY.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC
      USE GWFLPFMODULE,ONLY:HK,CHANI,HANI
C     ------------------------------------------------------------------
C
      ZERO=0.
      TWO=2.
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
C
C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HK(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
C
C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         T1=HK(J,I,K)*CC(J,I,K)
C3A-----IF THIS IS NOT THE LAST COLUMN (RIGHTMOST), CALCULATE
C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
               T2=HK(J+1,I,K)*CC(J+1,I,K)
               CR(J,I,K)=TWO*T2*T1*DELC(I)/(T1*DELR(J+1)+T2*DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
C3B-----IF THIS IS THE LAST COLUMN, SET BRANCH CONDUCTANCE=0.
            CR(J,I,K)=ZERO
         END IF
C
C3C-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
C3C-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
               T2=HK(J,I+1,K)*CC(J,I+1,K)
               IF(CHANI(K).LE.ZERO) THEN
                  KHANI=-CHANI(K)
                  T1=T1*HANI(J,I,KHANI)
                  T2=T2*HANI(J,I+1,KHANI)
               ELSE
                  T1=T1*CHANI(K)
                  T2=T2*CHANI(K)
               END IF
               CC(J,I,K)=TWO*T2*T1*DELR(J)/(T1*DELC(I+1)+T2*DELC(I))
            ELSE
C3D-----IF THIS IS THE LAST ROW, SET BRANCH CONDUCTANCE=0.
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
C
C4------RETURN
      RETURN
      END SUBROUTINE SGWF2LPF7HHARM

!***********************************************************************

      SUBROUTINE SGWF2LPF7HLOG(K)
C     ******************************************************************
C-----COMPUTE HORIZONTAL CONDUCTANCE USING LOGARITHMIC MEAN
C-----TRANSMISSIVITY -- ACTIVATED BY LAYAVG=1
C-----CELL SATURATED THICKNESS IS IN CC.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC
      USE GWFLPFMODULE,ONLY:HK,CHANI,HANI
C     ------------------------------------------------------------------
C
      ZERO=0.
      TWO=2.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
C
C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HK(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
C
C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         T1=HK(J,I,K)*CC(J,I,K)
C3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
C3A1----LOGARITHMIC MEAN INTERBLOCK TRANSMISSIVITY
               T2=HK(J+1,I,K)*CC(J+1,I,K)
               RATIO=T2/T1
               IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                  T=(T2-T1)/LOG(RATIO)
               ELSE
                  T=HALF*(T1+T2)
               END IF
               CR(J,I,K)=TWO*DELC(I)*T/(DELR(J+1)+DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
            CR(J,I,K)=ZERO
         END IF
C
C3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
C3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
               T2=HK(J,I+1,K)*CC(J,I+1,K)
               IF(CHANI(K).LE.ZERO) THEN
                  KHANI=-CHANI(K)
                  T1=T1*HANI(J,I,KHANI)
                  T2=T2*HANI(J,I+1,KHANI)
               ELSE
                  T1=T1*CHANI(K)
                  T2=T2*CHANI(K)
               END IF
               RATIO=T2/T1
               IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                  T=(T2-T1)/LOG(RATIO)
               ELSE
                  T=HALF*(T1+T2)
               END IF
               CC(J,I,K)=TWO*DELR(J)*T/(DELC(I+1)+DELC(I))
            ELSE
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
C
C4------RETURN
      RETURN
      END SUBROUTINE SGWF2LPF7HLOG

!***********************************************************************

      SUBROUTINE SGWF2LPF7HUNCNF(K)
C     ******************************************************************
C-----COMPUTE HORIZONTAL CONDUCTANCE USING ARITHMETIC MEAN SATURATED
C-----THICKNESS AND LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY.
C-----CELL SATURATED THICKNESS IS IN CC.
C-----ACTIVATED BY LAYAVG=2
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC
      USE GWFLPFMODULE,ONLY:HK,CHANI,HANI
C     ------------------------------------------------------------------
C
      ZERO=0.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
C
C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
C
C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HK(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
C
C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         HYC1=HK(J,I,K)
C3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
C3A1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
               HYC2=HK(J+1,I,K)
               RATIO=HYC2/HYC1
               IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                  HYC=(HYC2-HYC1)/LOG(RATIO)
               ELSE
                  HYC=HALF*(HYC1+HYC2)
               END IF
C3A2----MULTIPLY LOGARITHMIC K BY ARITMETIC SATURATED THICKNESS.
               CR(J,I,K)=DELC(I)*HYC*(CC(J,I,K)+CC(J+1,I,K))/
     1                      (DELR(J+1)+DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
            CR(J,I,K)=ZERO
         END IF
C
C3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
C3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
C3B1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
               HYC2=HK(J,I+1,K)
               IF(CHANI(K).LE.ZERO) THEN
                  KHANI=-CHANI(K)
                  HYC1=HYC1*HANI(J,I,KHANI)
                  HYC2=HYC2*HANI(J,I+1,KHANI)
               ELSE
                  HYC1=HYC1*CHANI(K)
                  HYC2=HYC2*CHANI(K)
               END IF
               RATIO=HYC2/HYC1
               IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                  HYC=(HYC2-HYC1)/LOG(RATIO)
               ELSE
                  HYC=HALF*(HYC1+HYC2)
               END IF
C3B2----MULTIPLY LOGARITHMIC K BY ARITMETIC SATURATED THICKNESS.
               CC(J,I,K)=DELR(J)*HYC*(CC(J,I,K)+CC(J,I+1,K))/
     1                      (DELC(I+1)+DELC(I))
            ELSE
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
C
C4------RETURN.
      RETURN
      END SUBROUTINE SGWF2LPF7HUNCNF

!***********************************************************************

      SUBROUTINE SGWF2LPF7VCOND(K)
C     ******************************************************************
C     COMPUTE VERTICAL BRANCH CONDUCTANCE BETWEEN A LAYER AND THE NEXT
C     LOWER LAYER FROM VERTICAL HYDRAULIC CONDUCTIVITY.
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,CV,DELR,DELC,
     1                        BOTM,LBOTM,LAYCBD,IOUT,STRT
!      USE GWFLPFMODULE,  ONLY:LAYTYP,LAYAVG,CHANI,LAYVKA,LAYWET,
!     1                        HK,VKA,VKCB,NOCVCO,ICONCV,LAYSTRT
      use GWFLPFMODULE, only: LAYTYP,LAYVKA,
     1                        HK,VKA,VKCB,NOCVCO,ICONCV,LAYSTRT
C
      DOUBLE PRECISION BBOT,TTOP,HHD
C     ------------------------------------------------------------------
C
      IF(K.EQ.NLAY) RETURN
      ZERO=0.
      HALF=0.5
C
C1------LOOP THROUGH ALL CELLS IN THE LAYER.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
      CV(J,I,K)=ZERO
      IF(IBOUND(J,I,K).NE.0 .AND. IBOUND(J,I,K+1).NE.0) THEN
C
C2------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL.
         IF(LAYVKA(K).EQ.0) THEN
            HYC1=VKA(J,I,K)
         ELSE
            HYC1=HK(J,I,K)/VKA(J,I,K)
         END IF
         IF(HYC1.GT.ZERO) THEN
C3------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL BELOW.
            IF(LAYVKA(K+1).EQ.0) THEN
               HYC2=VKA(J,I,K+1)
            ELSE
               HYC2=(HK(J,I,K+1)/VKA(J,I,K+1))
            END IF
            IF(HYC2.GT.ZERO) THEN
C
C4------CALCULATE INVERSE LEAKANCE FOR CELL.  ICONCV FLAG PREVENTS
C4------CV FROM BEING HEAD DEPENDENT.
               BBOT=BOTM(J,I,LBOTM(K))
               TTOP=BOTM(J,I,LBOTM(K)-1)
               IF(LAYSTRT(K).NE.0) TTOP=STRT(J,I,K)
               IF(LAYTYP(K).NE.0 .AND. ICONCV.EQ.0) THEN
                  HHD=HNEW(J,I,K)
                  IF(HHD.LT.TTOP) TTOP=HHD
               END IF
               BOVK1=(TTOP-BBOT)*HALF/HYC1
C
C5------CALCULATE INVERSE LEAKANCE FOR CELL BELOW.
               BBOT=BOTM(J,I,LBOTM(K+1))
               TTOP=BOTM(J,I,LBOTM(K+1)-1)
               IF(LAYSTRT(K+1).NE.0) TTOP=STRT(J,I,K+1)
               B=(TTOP-BBOT)*HALF
C
C5A-----IF CELL BELOW IS NOT SATURATED, DO NOT INCLUDE ITS CONDUCTANCE
C5A-----IN THE VERTICAL CONDUCTANCE CALULATION, EXCEPT THAT THE NOCVCO
C5A-----AND ICONCV FLAGS TURN OFF THIS CORRECTION.
               IF(LAYTYP(K+1).NE.0
     1                    .AND.NOCVCO.EQ.0 .AND. ICONCV.EQ.0) THEN
                  HHD=HNEW(J,I,K+1)
                  IF(HHD.LT.TTOP) B=ZERO
               END IF
               BOVK2=B/HYC2
C
C6------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CONFINING BED.
               IF(LAYCBD(K).NE.0) THEN
                  IF(VKCB(J,I,LAYCBD(K)).GT.ZERO) THEN
C
C7------CALCULATE INVERSE LEAKANCE FOR CONFINING BED.
                     B=BOTM(J,I,LBOTM(K))-BOTM(J,I,LBOTM(K)+1)
                     IF(B.LT.ZERO) THEN
                        WRITE(IOUT,45) K,I,J
   45                   FORMAT(1X,/1X,
     1  'Negative confining bed thickness below cell (Layer,row,col)',
     2                  I4,',',I5,',',I5)
            WRITE(IOUT,46) BOTM(J,I,LBOTM(K)),BOTM(J,I,LBOTM(K)+1)
   46       FORMAT(1X,'Top elevation, bottom elevation:',1P,2G13.5)
                        CALL USTOP(' ')
                     END IF
                     CBBOVK=B/VKCB(J,I,LAYCBD(K))
                     CV(J,I,K)=DELR(J)*DELC(I)/(BOVK1+CBBOVK+BOVK2)
                  END IF
               ELSE
                  CV(J,I,K)=DELR(J)*DELC(I)/(BOVK1+BOVK2)
               END IF
            END IF
         END IF
      END IF
  100 CONTINUE
C
C8------RETURN.
      RETURN
      END SUBROUTINE SGWF2LPF7VCOND

!***********************************************************************

      SUBROUTINE SGWF2LPF7CK(IOUT,NP,PTYP)
C     ******************************************************************
C     CHECK THAT JUST-DEFINED PARAMETER OF TYPE 'VK' OR 'VANI' IS USED
C     CONSISTENTLY WITH LAYVKA ENTRIES FOR LAYERS LISTED IN CLUSTERS FOR
C     THE PARAMETER
C     ******************************************************************
C
C      SPECIFICATIONS:
C     ------------------------------------------------------------------
      !USE GWFLPFMODULE,  ONLY:LAYTYP,LAYAVG,CHANI,LAYVKA,LAYWET
      use GWFLPFMODULE, only: layvka
      USE PARAMMODULE
C
      CHARACTER*4 PTYP
C     ------------------------------------------------------------------
C
C1------LOOP THROUGH THE CLUSTERS FOR THIS PARAMETER.
      DO 10 ICL = IPLOC(1,NP),IPLOC(2,NP)
        LAY = IPCLST(1,ICL)
        LV = LAYVKA(LAY)
        IF (PTYP.EQ.'VK  ' .AND. LV.NE.0) THEN
          WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VK'
  590     FORMAT(/,
     &1X,'LAYVKA entered for layer ',i3,' is: ',i3,'; however,',
     &' layer ',i3,' is',/,' listed in a cluster for parameter "',a,
     &'" of type ',a,' and')
          WRITE (IOUT,600)
  600     FORMAT(
     &1X,'parameters of type VK can apply only to layers for which',
     &/,' LAYVKA is specified as zero -- STOP EXECUTION (SGWF2LPF7CK)')
          CALL USTOP(' ')
        ELSEIF (PTYP.EQ.'VANI' .AND. LV.EQ.0) THEN
          WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VANI'
          WRITE (IOUT,610)
  610     FORMAT(
     &1X,'parameters of type VANI can apply only to layers for which',/,
     &' LAYVKA is not specified as zero -- STOP EXECUTION',
     &' (SGWF2LPF7CK)')
          CALL USTOP(' ')
        ENDIF
   10 CONTINUE
C
C2------Return.
      RETURN
      END SUBROUTINE SGWF2LPF7CK

      end module GwfLpfSubs

