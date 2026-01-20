module UpwSubsModule
  use ConstantsPHMFModule, only: HDRYDEFAULT
  use ModelModule, only: ModelType
  use SimPHMFModule, only: store_error, store_note, store_warning

  private
  public :: GWF2UPW1AR

contains

      SUBROUTINE GWF2UPW1AR(In, model)

      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD, &
                          NCNFBD,IBOUND,BUFF,BOTM,NBOTM,DELR,DELC,IOUT, &
                          LBOTM,HNEW
      USE GWFBASMODULE,ONLY:HDRY
      USE GWFNWTMODULE, ONLY: Numcell
      USE GWFUPWMODULE
      use SimPHMFModule, only: ustop
      use utl7module, only: URDCOM, URWORD, U2DREL
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
      INTRINSIC INT
      !EXTERNAL URDCOM, URWORD
      !EXTERNAL SGWF2UPW1PSV
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER In
      type(ModelType), pointer, intent(inout) :: model
!     ------------------------------------------------------------------
!     LOCAL VARIABLES
!     ------------------------------------------------------------------
      INTEGER lloc, istart, istop, i, ic, ir, il, jj
      CHARACTER(LEN=200) line
      INTEGER NPHK,NPVKCB,NPVK,NPVANI,NPSS,NPSY,NPHANI
      INTEGER IANAME,KHANI,N,KK,j,k,NCNVRT,NHANI,NWETD
!     LOCAL VARIABLES FOR DEFINING CELL PROPERTES (FROM LPF)
      INTEGER NPUPW, NOPCHK
      double precision :: ZERO, R
!
      CHARACTER*14 LAYPRN(5),AVGNAM(3),TYPNAM(2),VKANAM(2),WETNAM(2), &
                 HANNAM
      DATA AVGNAM/'      HARMONIC','   LOGARITHMIC','     LOG-ARITH'/
      DATA TYPNAM/'      CONFINED','   CONVERTIBLE'/
      DATA VKANAM/'    VERTICAL K','    ANISOTROPY'/
      DATA WETNAM/'  NON-WETTABLE','      WETTABLE'/
      DATA HANNAM/'      VARIABLE'/
      CHARACTER*24 ANAME(9),STOTXT
      CHARACTER*4 PTYP
      character(len=1000) :: msg
      integer :: ianyharmonic, ianylog, ianyarithlog, iavg, knew
      logical :: NovfcFound, ConstantCvFound, NoCvCorrectionFound
      logical :: Note10Written = .false.
!C
      DATA ANAME(1) /'   HYD. COND. ALONG ROWS'/
      DATA ANAME(2) /'  HORIZ. ANI. (COL./ROW)'/
      DATA ANAME(3) /'     VERTICAL HYD. COND.'/
      DATA ANAME(4) /' HORIZ. TO VERTICAL ANI.'/
      DATA ANAME(5) /'QUASI3D VERT. HYD. COND.'/
      DATA ANAME(6) /'        SPECIFIC STORAGE'/
      DATA ANAME(7) /'          SPECIFIC YIELD'/
      DATA ANAME(9) /'     STORAGE COEFFICIENT'/
      ! formats
   10 format('In UPW, storage changes under unconfined conditions ', &
      'account for both specific yield and specific storage.', &
      ' As a consequence, numerical results of MODFLOW 6 may', &
      ' differ from the results of MODFLOW-NWT.')
!     ------------------------------------------------------------------
!1------IDENTIFY PACKAGE AND INITIALIZE.
      WRITE (Iout, 9001) In
 9001 FORMAT (1X, /' UPW1 -- UPSTREAM WEIGHTING FLOW PACKAGE, ', &
            'VERSION 1.0.6, 12/05/2012', /, 9X, 'INPUT READ FROM UNIT', &
             I3,/)
!  ALLOCATE, READ AND SET DATA FOR CELL PROPERTIES (FROM LPF)
!C1------Allocate scalar data.
      ALLOCATE(IUPWCB,NOVFC,Iuupw)
!  STORE UPW UNIT NUMBER IN MODULE VARIABLE
      Iuupw = In
      ALLOCATE(ISFAC,ICONCV,ITHFLG,NOCVCO,IPHDRY)
      ZERO=0.
      write(*,*)'Processing UPW package input...'
      call model%NpfWriter%AllocateArrays()
      model%Newton = .true.
!C
!C3------READ COMMENTS AND ITEM 1.
      CALL URDCOM(IN,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IUPWCB,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,HDRY,IOUT,IN)
      hdry = hdrydefault
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,NPUPW,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IPHDRY,R,IOUT,IN)
      model%NpfWriter%Inpfcb = IUPWCB
      model%NpfWriter%Hdry = HDRY
!
!C
!C3A-----WRITE ITEM 1
      IF(IUPWCB.LT.0) WRITE(IOUT,8)
    8 FORMAT(1X,'CONSTANT-HEAD CELL-BY-CELL FLOWS WILL BE PRINTED', &
       ' WHEN ICBCFL IS NOT 0')
      IF(IUPWCB.GT.0) WRITE(IOUT,9) IUPWCB
    9 FORMAT(1X,'CELL-BY-CELL FLOWS WILL BE SAVED ON UNIT ',I4)
      IF(NPUPW.GT.0) THEN
         WRITE(IOUT,15) NPUPW
   15    FORMAT(1X,I5,' Named Parameters     ')
      ELSE
         NPUPW=0
         WRITE(IOUT,'(A)') ' No named parameters'
      END IF
!C
!C3B-----GET OPTIONS.
      ISFAC=0
      ICONCV=0
      ITHFLG=0
      NOCVCO=0
      NOVFC=0
      NOPCHK=0
      STOTXT=ANAME(6)
      !
      ! -- Record UPW options as specified
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
   21    FORMAT(1X,'STORAGECOEFFICIENT OPTION:',/, &
          1X,'Read storage coefficient rather than specific storage' &
          1X,'Option not supported in UPW Package for MODFLOW-NWT')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'CONSTANTCV') THEN
         ConstantCvFound = .true.
         ICONCV=1
         WRITE(IOUT,23)
   23    FORMAT(1X,'CONSTANTCV OPTION:',/,1X,'Constant vertical', &
              ' conductance for convertible layers' &
          1X,'Option not supported in UPW Package for MODFLOW-NWT')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'THICKSTRT') THEN
         ITHFLG=1
         WRITE(IOUT,25)
   25    FORMAT(1X,'THICKSTRT OPTION:',/,1X,'Negative LAYTYP indicates', &
      ' confined layer with thickness computed from STRT-BOT' &
          1X,'Option not supported in UPW Package for MODFLOW-NWT')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOCVCORRECTION') THEN
         NoCvCorrectionFound = .true.
         NOCVCO=1
         WRITE(IOUT,27)
   27    FORMAT(1X,'NOCVCORRECTION OPTION:',/,1X, &
         'Do not adjust vertical conductance when applying', &
                   ' the vertical flow correction' &
          1X,'Option not supported, no vertical flow correction ' &
             ' in UPW Package')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOVFC') THEN
         NovfcFound = .true.
         NOVFC=1
         NOCVCO=1
         WRITE(IOUT,29)
   29    FORMAT(1X,'NOVFC OPTION:',/,1X, &
         'vertical flow correction does not apply in UPW Package')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'NOPARCHECK') THEN
         NOPCHK=1
         WRITE(IOUT,30)
   30    FORMAT(1X,'NOPARCHECK  OPTION:',/,1X, &
         'For data defined by parameters, do not check to see if ', &
             'parameters define data at all cells')
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
!C
!C4------ALLOCATE AND READ LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET, LAYSTRT.
      ALLOCATE (Sn(Numcell), So(Numcell))
      Sn = 0.0D0
      So = 0.0D0
      ALLOCATE(LAYTYPUPW(NLAY))
      ALLOCATE(LAYAVG(NLAY))
      ALLOCATE(CHANI(NLAY))
      ALLOCATE(LAYVKAUPW(NLAY))
      ALLOCATE(LAYWET(NLAY))
      ALLOCATE(LAYSTRT(NLAY))
      ALLOCATE(IBOUND2(NCOL,NROW,NLAY))
      IBOUND2 = 0
      READ(IN,*) (LAYTYPUPW(K),K=1,NLAY)
      READ(IN,*) (LAYAVG(K),K=1,NLAY)
      READ(IN,*) (CHANI(K),K=1,NLAY)
      READ(IN,*) (LAYVKAUPW(K),K=1,NLAY)
      READ(IN,*) (LAYWET(K),K=1,NLAY)
      ! Note: According to UPW input instructions, LAYWET
      !       should always be zero.
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
        select case (laytypUPW(k))
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
          model%NpfWriter%NumConvertible = &
                 model%NpfWriter%NumConvertible + 1
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
          msg = ' Non-uniform interblock-transmissivity averaging ' &
              // 'methods found in LPF input. HARMONIC is used in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'HARMONIC'
        elseif (ianylog == 1) then
          msg = ' Non-uniform interblock-transmissivity averaging ' &
              // 'methods found in LPF input. LOGARITHMIC is used ' &
              // 'in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'LOGARITHMIC'
        elseif (ianyarithlog == 1) then
          msg = ' Non-uniform interblock-transmissivity averaging ' &
              // 'methods found in LPF input. AMT-LMK is used ' &
              // 'in NPF.'
          call store_warning(msg)
          model%NpfWriter%CellAveraging = 'AMT-LMK'
        endif
      endif
      !
!C
!C4A-----PRINT A TABLE OF VALUES FOR LAYTYP, LAYAVG, CHANI, LAYVKA, LAYWET.
      WRITE(IOUT,47)
   47 FORMAT(1X,/3X,'LAYER FLAGS:',/1X, &
      'LAYER       LAYTYP          LAYAVG    CHANI    ', &
      '       LAYVKA           LAYWET',/1X,75('-'))
      DO 50 K=1,NLAY
      WRITE(IOUT,48) K,LAYTYPUPW(K),LAYAVG(K),CHANI(K),LAYVKAUPW(K), &
                    LAYWET(K)
   48 FORMAT(1X,I4,2I14,1PE14.3,2I14)
!C
!C4A1----SET GLOBAL HEAD-DEPENDENT TRANSMISSIVITY AND STORAGE FLAGS.
      IF (LAYTYPUPW(K).GT.0) THEN
        LAYHDT(K)=1
        LAYHDS(K)=1
      ELSE
        LAYHDT(K)=0
        LAYHDS(K)=0
      ENDIF
   50 CONTINUE
!C
!C4A2----SET LAYSTRT AND RESET LAYTYP IF THICKSTRT OPTION IS ACTIVE.
      DO 60 K=1,NLAY
      LAYSTRT(K)=0
      IF(LAYTYPUPW(K).LT.0 .AND. ITHFLG.NE.0) THEN
         LAYSTRT(K)=1
         LAYTYPUPW(K)=0
         LAYHDT(K)=0
         LAYHDS(K)=0
         WRITE(IOUT,57) K
   57    FORMAT(1X,'Layer',I5, &
      ' is confined because LAYTYP<0 and THICKSTRT option is active')
      END IF
   60 CONTINUE
!C
!C4B-----BASED ON LAYTYP, LAYAVG, CHANI, LAYWET, COUNT THE NUMBER OF EACH
!C4B-----TYPE OF 2-D ARRAY; CHECK VALUES FOR CONSISTENCY; AND SETUP
!C4B-----POINTERS IN LAYTYP, CHANI, AND LAYWET FOR CONVENIENT ACCESS
!C4B-----TO SC2, HANI, and WETDRY.  PRINT INTERPRETED VALUES OF FLAGS.
      NCNVRT=0
      NHANI=0
      NWETD=0
      WRITE(IOUT,67)
   67 FORMAT(1X,/3X,'INTERPRETATION OF LAYER FLAGS:',/1X, &
       '                       INTERBLOCK     HORIZONTAL', &
       '    DATA IN',/1X, &
       '        LAYER TYPE   TRANSMISSIVITY   ANISOTROPY', &
       '   ARRAY VKA   WETTABILITY',/1X, &
       'LAYER   (LAYTYP)        (LAYAVG)      (CHANI)   ', &
       '  (LAYVKA)       (LAYWET)',/1X,75('-'))
      DO 100 K=1,NLAY
      IF(LAYTYPUPW(K).GT.0) THEN
         NCNVRT=NCNVRT+1
         LAYTYPUPW(K)=NCNVRT
      END IF
      IF(CHANI(K).LE.ZERO) THEN
         NHANI=NHANI+1
         CHANI(K)=-NHANI
      END IF
      IF(LAYWET(K).NE.0) THEN
!         IF(LAYTYPUPW(K).EQ.0) THEN
            WRITE(IOUT,*)
            WRITE(IOUT,*) &
               ' LAYWET is not 0 and wetting does not apply in UPW '
            WRITE(IOUT,*) ' LAYWET must be 0 when using the UPW Package'
            CALL USTOP(' ')
!         ELSE
!            NWETD=NWETD+1
!            LAYWET(K)=NWETD
!         END IF
      END IF
      IF(LAYAVG(K).LT.0 .OR. LAYAVG(K).GT.2) THEN
         WRITE(IOUT,74) LAYAVG(K)
   74    FORMAT(1X,I8, &
         ' IS AN INVALID LAYAVG VALUE -- MUST BE 0, 1, or 2')
         CALL USTOP(' ')
      END IF
      LAYPRN(1)=TYPNAM(1)
      IF(LAYTYPUPW(K).GT.0) LAYPRN(1)=TYPNAM(2)
      LAYPRN(2)=AVGNAM(LAYAVG(K)+1)
      IF(CHANI(K).LE.0) THEN
         LAYPRN(3)=HANNAM
      ELSE
         WRITE(LAYPRN(3),'(1PE14.3)') CHANI(K)
      END IF
      LAYPRN(4)=VKANAM(1)
      IF(LAYVKAUPW(K).NE.0) LAYPRN(4)=VKANAM(2)
      LAYPRN(5)=WETNAM(1)
      IF(LAYWET(K).NE.0) LAYPRN(5)=WETNAM(2)
      WRITE(IOUT,78) K,(LAYPRN(I),I=1,5)
   78 FORMAT(1X,I4,5A)
  100 CONTINUE
!C
!C4C-----PRINT WETTING INFORMATION.  RGN commented out because this does not apply
!      IF(NWETD.EQ.0) THEN
!         WRITE(IOUT,13)
!   13    FORMAT(1X,/,1X,'WETTING CAPABILITY IS NOT ACTIVE IN ANY LAYER')
!         IWDFLG=0
!      ELSE
!         WRITE(IOUT,12) NWETD
!   12    FORMAT(1X,/,1X,'WETTING CAPABILITY IS ACTIVE IN',I4,' LAYERS')
!         IWDFLG=1
!         READ(IN,*) WETFCT,IWETIT,IHDWET
!         IF(IWETIT.LE.0) IWETIT=1
!         WRITE(IOUT,*) ' WETTING FACTOR=',WETFCT
!         WRITE(IOUT,*) ' WETTING ITERATION INTERVAL=',IWETIT
!         WRITE(IOUT,*) ' IHDWET=',IHDWET
!      END IF
! ALLOCATE THESE BECAUSE THEY ARE LEFT IN THE CODE (DEALLOCATED)
      ALLOCATE(WETFCT,IWETIT,IHDWET,IWDFLG)
      ! Not sure if it matters if REWET is listed as an option or not.
      ! I assume REWET is redundant when UPW is used.
      model%NpfWriter%Rewet = .false.
!C
!C5------ALLOCATE MEMORY FOR ARRAYS.
      ALLOCATE(LAYFLG(6,NLAY))
      ALLOCATE(HKUPW(NCOL,NROW,NLAY))
      ALLOCATE(VKAUPW(NCOL,NROW,NLAY))
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
      SC1 = 0.0
! RGN 6/25/09      IF(ITRSS.NE.0 .AND. NCNVRT.GT.0) THEN
      IF(ITRSS.NE.0) THEN
! RGN 6/25/09        ALLOCATE(SC2UPW(NCOL,NROW,NCNVRT))
         ALLOCATE(SC2UPW(NCOL,NROW,NLAY))
      ELSE
         ALLOCATE(SC2UPW(1,1,1))
      END IF
      SC2UPW = 0.0
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
!C
!C6------READ PARAMETER DEFINITIONS
      NPHK=0
      NPVKCB=0
      NPVK=0
      NPVANI=0
      NPSS=0
      NPSY=0
      NPHANI=0
      IF(NPUPW.GT.0) THEN
         WRITE(IOUT,115)
  115    FORMAT(/,' PARAMETERS DEFINED IN THE UPW PACKAGE')
         DO 120 K=1,NPUPW
         CALL UPARARRRP(IN,IOUT,N,1,PTYP,1,0,-1)
!C   Note that NPHK and the other NP variables in
!C   this group are used only as flags, not counts
         IF(PTYP.EQ.'HK') THEN
            NPHK=1
         ELSE IF(PTYP.EQ.'HANI') THEN
!C6A-----WHEN A HANI PARAMETER IS USED, THEN ALL HORIZONTAL ANISOTROPY
!C6A-----MUST BE DEFINED USING PARAMETERS.  ENSURE THAT ALL CHANI <= 0
            DO 118 I = 1, NLAY
              IF (CHANI(I).GT.0.0) THEN
                WRITE(IOUT,117)
  117           FORMAT(/, &
     'ERROR: WHEN A HANI PARAMETER IS USED, CHANI FOR ALL LAYERS',/, &
     'MUST BE LESS THAN OR EQUAL TO 0.0 -- STOP EXECUTION', &
     ' (GWF2UPW1AR)')
                CALL USTOP(' ')
              ENDIF
  118       CONTINUE
            NPHANI=1
         ELSE IF(PTYP.EQ.'VKCB') THEN
            NPVKCB=1
         ELSE IF(PTYP.EQ.'VK') THEN
            NPVK=1
            CALL SGWF2UPWCK(IOUT,N,'VK  ')
         ELSE IF(PTYP.EQ.'VANI') THEN
            NPVANI=1
            CALL SGWF2UPWCK(IOUT,N,'VANI')
         ELSE IF(PTYP.EQ.'SS') THEN
            NPSS=1
         ELSE IF(PTYP.EQ.'SY') THEN
            NPSY=1
         ELSE
            WRITE(IOUT,*) ' Invalid parameter type for UPW Package'
            CALL USTOP(' ')
         END IF
  120    CONTINUE
      END IF
!C
!C7------DEFINE DATA FOR EACH LAYER -- VIA READING OR NAMED PARAMETERS.
      knew = 0
      DO 200 K=1,NLAY
      knew = knew + 1
      KK=K
!C
!C7A-----DEFINE HORIZONTAL HYDRAULIC CONDUCTIVITY (HK)
      IF(NPHK.EQ.0) THEN
         CALL U2DREL(HKUPW(:,:,KK),ANAME(1),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(1,K)
         WRITE(IOUT,121) ANAME(1),K,LAYFLG(1,K)
  121    FORMAT(1X,/1X,A,' FOR LAYER',I4, &
        ' WILL BE DEFINED BY PARAMETERS',/1X,'(PRINT FLAG=',I4,')')
         CALL UPARARRSUB1(HKUPW(:,:,KK),NCOL,NROW,KK,'HK', &
           IOUT,ANAME(1),LAYFLG(1,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,NLAY, &
           NROW,'HK  ')
      END IF
      !
      ! Copy HK to NpfWriter%HK -- This accounts for parameter values
      ! read from PVAL file, and zone and mult arrays.
      model%NpfWriter%HK(:,:,K) = HKUPW(:,:,KK)
      !
!C
!C7B-----READ HORIZONTAL ANISOTROPY IF CHANI IS NON-ZERO
      IF(CHANI(K).LE.ZERO) THEN
        KHANI=-CHANI(K)
        IF(NPHANI.EQ.0) THEN
           CALL U2DREL(HANI(:,:,KHANI),ANAME(2),NROW,NCOL,KK,IN,IOUT)
        ELSE
           READ(IN,*) LAYFLG(6,K)
           WRITE(IOUT,121) ANAME(2),K,LAYFLG(6,K)
           CALL UPARARRSUB1(HANI(:,:,KHANI),NCOL,NROW,KK,'HANI', &
           IOUT,ANAME(2),LAYFLG(6,KK))
           IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL, &
           NLAY,NROW,'HANI')
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
!C
!C7C-----DEFINE VERTICAL HYDRAULIC CONDUCTIVITY OR HORIZONTAL TO VERTICAL
!C7C-----ANISOTROPY (VKA).
      IANAME=3
      PTYP='VK'
      IF(LAYVKAUPW(K).NE.0) THEN
         IANAME=4
         PTYP='VANI'
      END IF
      IF(NPVK.EQ.0 .AND. NPVANI.EQ.0) THEN
         CALL U2DREL(VKAUPW(:,:,KK),ANAME(IANAME),NROW,NCOL,KK,IN,IOUT)
      ELSE
         READ(IN,*) LAYFLG(2,K)
         WRITE(IOUT,121) ANAME(IANAME),K,LAYFLG(2,K)
         CALL UPARARRSUB1(VKAUPW(:,:,KK),NCOL,NROW,KK,PTYP,IOUT, &
                            ANAME(IANAME),LAYFLG(2,KK))
         IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL,NLAY, &
                            NROW,PTYP)
      END IF
      !
      ! Define VK for mf6
      if (LAYVKAUPW(k) == 0) then
        ! VKA contains vertical hydraulic conductivity
        model%NpfWriter%vk(:,:,knew) = VKAUPW(:,:,k)
      else
        ! VKA contains vertical anisotropy. Compute VK from HK and anisotropy.
        ! Anisotropy is ratio HK/VK.
        do i=1,NROW
          do j=1,NCOL
            model%NpfWriter%vk(j,i,knew) = HKUPW(j,i,k) / VKAUPW(j,i,k)
          enddo
        enddo
      endif
!C
!C7D-----DEFINE SPECIFIC STORAGE OR STORAGE COEFFICIENT IN ARRAY SC1 IF TRANSIENT.
      IF(ITRSS.NE.0) THEN
         IF(NPSS.EQ.0) THEN
            CALL U2DREL(SC1(:,:,KK),STOTXT,NROW,NCOL,KK,IN,IOUT)
         ELSE
            READ(IN,*) LAYFLG(3,K)
            WRITE(IOUT,121) STOTXT,K,LAYFLG(3,K)
            CALL UPARARRSUB1(SC1(:,:,KK),NCOL,NROW,KK,'SS', &
                IOUT,STOTXT,LAYFLG(3,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL, &
                NLAY,NROW,'SS  ')
         END IF
         ! For current layer, define Ss array of NpfWriter
         ! here (before Ss is multiplied by cell thickness
         ! and SC1 is multiplied by cell area).
         model%StoWriter%Ss(:,:,k) = sc1(:,:,k)
         IF(ISFAC.EQ.0) THEN
            CALL SGWF2UPWSC(SC1(:,:,KK),KK,1)
         ELSE
            CALL SGWF2UPWSC(SC1(:,:,KK),KK,0)
         END IF
      END IF
!C
!C7E-----DEFINE SPECIFIC YIELD IN ARRAY SC2 IF TRANSIENT AND LAYER IS
!C7E-----IS CONVERTIBLE.
      IF(LAYTYPUPW(K).GT.0) THEN
         IF(ITRSS.NE.0) THEN
            IF(NPSY.EQ.0) THEN
             CALL U2DREL(SC2UPW(:,:,LAYTYPUPW(K)),ANAME(7),NROW,NCOL,KK, &
                      IN,IOUT)
            ELSE
               READ(IN,*) LAYFLG(4,K)
               WRITE(IOUT,121) ANAME(7),K,LAYFLG(4,K)
               CALL UPARARRSUB1(SC2UPW(:,:,LAYTYPUPW(K)),NCOL, &
                   NROW,KK,'SY',IOUT,ANAME(7),LAYFLG(4,KK))
               IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K, &
                   NCOL,NLAY,NROW,'SY  ')
            END IF
            ! For current layer, define Sy array of NpfWriter
            ! here (before Sy is multiplied by cell area).
            model%StoWriter%Sy(:,:,k) = SC2UPW(:,:,LAYTYPUPW(k))
            CALL SGWF2UPWSC(SC2UPW(:,:,LAYTYPUPW(K)),KK,0)
         END IF
      ELSE
        IF(ITRSS.NE.0) THEN
          DO J=1,NROW
            DO I=1,NCOL
              SC2UPW(I,J,KK) = 0.0D0    !SC1(I,J,KK)
            END DO
          END DO
        END IF
      END IF
!C
!C7F-----READ CONFINING BED VERTICAL HYDRAULIC CONDUCTIVITY (VKCB)
      IF(LAYCBD(K).NE.0) THEN
         IF(NPVKCB.EQ.0) THEN
            CALL U2DREL(VKCB(:,:,LAYCBD(K)),ANAME(5),NROW,NCOL,KK,IN, &
                  IOUT)
         ELSE
            READ(IN,*) LAYFLG(5,K)
            WRITE(IOUT,121) ANAME(5),K,LAYFLG(5,K)
            CALL UPARARRSUB1(VKCB(:,:,LAYCBD(K)),NCOL,NROW,KK, &
              'VKCB',IOUT,ANAME(5),LAYFLG(5,KK))
            IF(NOPCHK.EQ.0) CALL UPARARRCK(BUFF,IBOUND,IOUT,K,NCOL, &
              NLAY,NROW,'VKCB')
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
!C
!C7G-----READ WETDRY CODES IF WETTING CAPABILITY HAS BEEN INVOKED
!C7G-----(LAYWET NOT 0).
      IF(LAYWET(K).NE.0) THEN
         CALL U2DREL(WETDRY(:,:,LAYWET(K)),ANAME(8),NROW,NCOL,KK,IN, &
                 IOUT)
      END IF
  200 CONTINUE
!C
!C8------PREPARE AND CHECK LPF DATA.
      CALL SGWF2UPWN()
!C9------Calculate constant part of conductance. Conductance includes
!C       cell thickness for confined conditions.
      DO K=1,NLAY
        IF(LAYAVG(K).EQ.0) THEN
          IF ( LAYTYPUPW(K).GT.0 ) THEN
            CALL SGWF2UPW1HHARM(K)
          ELSE
            CALL SGWF2UPW1HHARMCON(K)
          END IF
        ELSE IF(LAYAVG(K).EQ.1) THEN
          IF ( LAYTYPUPW(K).GT.0 ) THEN
            CALL SGWF2UPW1HLOG(K)
          ELSE
            CALL SGWF2UPW1HLOGCON(K)
          END IF
        ELSE IF(LAYAVG(K).EQ.2) THEN
          IF ( LAYTYPUPW(K).GT.0 ) THEN
            CALL SGWF2UPW1HUNCNF(K)
          ELSE
            CALL SGWF2UPW1HUNCNFCON(K)
          END IF
        END IF
        CALL SGWF2UPW1VCOND(K)
      END DO
!10-----SAVE POINTERS FOR GRID AND RETURN
      CALL SGWF2UPW1PSV(model%IGrid)
!
!11-----RETURN
      END SUBROUTINE GWF2UPW1AR

      SUBROUTINE SGWF2UPWN()
!C     ******************************************************************
!C     INITIALIZE AND CHECK UPW DATA
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,LAYCBD,CV, &
                            BOTM,NBOTM,DELR,DELC,IOUT
      USE GWFBASMODULE,ONLY:HNOFLO
      USE GWFUPWMODULE,ONLY:HKUPW,VKAUPW
!C     ------------------------------------------------------------------
!C
!C1------DEFINE CONSTANTS.
      ZERO=0.
      HCNV=HNOFLO
!C
!C2------INSURE THAT EACH ACTIVE CELL HAS AT LEAST ONE NON-ZERO
!C2------TRANSMISSIVE PARAMETER.
      DO 60 K=1,NLAY
         DO 40 I=1,NROW
         DO 40 J=1,NCOL
         IF(IBOUND(J,I,K).EQ.0 ) GO TO 40
!C
!C3A-----CHECK HORIZONTAL HYDRAULIC CONDUCTIVITY (HK).
         IF(HKUPW(J,I,K).NE.ZERO) GO TO 40
!C
!C3B-----CHECK VERTICAL HYDRAULIC CONDUCTIVITY AND CONFINING BED
!C3B-----VERTICAL HYDRAULIC CONDUCTIVITY.
         IF(NLAY.GT.1) THEN

            IF(VKAUPW(J,I,K).NE.ZERO) THEN
               IF(K.NE.NLAY) THEN
                  IF (VKAUPW(J,I,K+1).NE.ZERO) GO TO 40
               END IF
               IF(K.NE.1) THEN
                  IF (VKAUPW(J,I,K-1).NE.ZERO) GO TO 40
               END IF
            END IF
         END IF
!C
!C3C-----ALL TRANSMISSIVE TERMS ARE ALL 0, SO CONVERT CELL TO NO FLOW.
         IBOUND(J,I,K)=0
         HNEW(J,I,K)=HCNV
         WRITE(IOUT,43) K,I,J
   40    CONTINUE
   43    FORMAT(1X,'NODE (LAYER,ROW,COL) ',I3,2(1X,I5), &
       ' ELIMINATED BECAUSE ALL HYDRAULIC',/, &
       ' CONDUCTIVITIES TO NODE ARE 0')
!C
   60 CONTINUE
!C
!C7------RETURN.
      RETURN
      END SUBROUTINE SGWF2UPWN

      SUBROUTINE SGWF2UPWSC(SC,K,ISPST)
!C     ******************************************************************
!C     COMPUTE STORAGE CAPACITY
!C     ******************************************************************
!C
!C     SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,DELR,DELC,BOTM,LBOTM,LAYCBD
!C
      double precision SC(NCOL,NROW)
!C     ------------------------------------------------------------------
!C
!C1------MULTIPLY SPECIFIC STORAGE BY THICKNESS, DELR, AND DELC TO GET
!C1------CONFINED STORAGE CAPACITY.
      IF(ISPST.NE.0) THEN
         DO 80 I=1,NROW
         DO 80 J=1,NCOL
! RGN Made this consistent with unconfined storage by not multiplying by thickness.
!         THICK=BOTM(J,I,LBOTM(K)-1)-BOTM(J,I,LBOTM(K))
!         SC(J,I)=SC(J,I)*THICK*DELR(J)*DELC(I)
         SC(J,I)=SC(J,I)*DELR(J)*DELC(I)
   80    CONTINUE
      ELSE
!C
!C2------MULTIPLY SPECIFIC YIELD BY DELR AND DELC TO GET UNCONFINED
!C2------STORAGEE CAPACITY(SC2).
         DO 85 I=1,NROW
         DO 85 J=1,NCOL
         SC(J,I)=SC(J,I)*DELR(J)*DELC(I)
   85    CONTINUE
      END IF
!C
      RETURN
      END SUBROUTINE SGWF2UPWSC

      SUBROUTINE SGWF2UPWCK(IOUT,NP,PTYP)
!C     ******************************************************************
!C     CHECK THAT JUST-DEFINED PARAMETER OF TYPE 'VK' OR 'VANI' IS USED
!C     CONSISTENTLY WITH LAYVKA ENTRIES FOR LAYERS LISTED IN CLUSTERS FOR
!C     THE PARAMETER
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GWFUPWMODULE,  ONLY:CHANI,LAYVKAUPW
      USE PARAMMODULE
      use SimPHMFModule, only: ustop
!C
      CHARACTER*4 PTYP
!C     ------------------------------------------------------------------
!C
!C1------LOOP THROUGH THE CLUSTERS FOR THIS PARAMETER.
      DO 10 ICL = IPLOC(1,NP),IPLOC(2,NP)
        LAY = IPCLST(1,ICL)
        LV = LAYVKAUPW(LAY)
        IF (PTYP.EQ.'VK  ' .AND. LV.NE.0) THEN
          WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VK'
  590     FORMAT(/, &
     1X,'LAYVKA entered for layer ',i3,' is: ',i3,'; however,', &
     ' layer ',i3,' is',/,' listed in a cluster for parameter "',a, &
     '" of type ',a,' and')
          WRITE (IOUT,600)
  600     FORMAT( &
     1X,'parameters of type VK can apply only to layers for which', &
     /,' LAYVKA is specified as zero -- ', &
       'STOP EXECUTION (SGWF2UPWCK)')
          CALL USTOP(' ')
        ELSEIF (PTYP.EQ.'VANI' .AND. LV.EQ.0) THEN
          WRITE (IOUT,590) LAY,LV,LAY,PARNAM(NP),'VANI'
          WRITE (IOUT,610)
  610     FORMAT( &
     1X,'parameters of type VANI can apply only to layers for which',/, &
     ' LAYVKA is not specified as zero -- STOP EXECUTION', &
     ' (SGWF2UPWCK)')
          CALL USTOP(' ')
        ENDIF
   10 CONTINUE
!C
!C2------Return.
      RETURN
      END SUBROUTINE SGWF2UPWCK

      SUBROUTINE SGWF2UPW1HHARM(K)
!C     ******************************************************************
!C      COMPUTE THE CONSTANT PART OF HORIZONTAL CONDUCTANCE BASED ON THE
!C      HARMONIC AVEARGE K FOR ADJACENT CELLS.
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, iout, &
                      CC, CR
      USE GWFUPWMODULE
!C     ------------------------------------------------------------------
!C
      ZERO=0.
      TWO=2.
!C
!C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
!C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
!C
!C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
!C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
         T1=HKUPW(J,I,K)
!C3A-----IF THIS IS NOT THE LAST COLUMN (RIGHTMOST), CALCULATE
!C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
               T2=HKUPW(J+1,I,K)
               CR(J,I,K)=TWO*T2*T1*DELC(I)/(T1*DELR(J+1)+T2*DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
!C3B-----IF THIS IS THE LAST COLUMN, SET BRANCH CONDUCTANCE=0.
            CR(J,I,K)=ZERO
         END IF
!C
!C3C-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
!C3C-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
               T2=HKUPW(J,I+1,K)
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
!C3D-----IF THIS IS THE LAST ROW, SET BRANCH CONDUCTANCE=0.
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
!C
!C4------RETURN
      RETURN
      END SUBROUTINE SGWF2UPW1HHARM

      SUBROUTINE SGWF2UPW1HLOG(K)
!C     ******************************************************************
!C-----COMPUTE CONSTANT PART OF HORIZONTAL CONDUCTANCE USING LOGARITHMIC
!C-----MEAN HYDRAULIC CONDUCTIVITY -- ACTIVATED BY LAYAVG=1
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, CR, CC
      USE GWFUPWMODULE
!C     ------------------------------------------------------------------
!C
      ZERO=0.
      TWO=2.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
!C
!C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
!C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
!C
!C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
!C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
!C
         T1=HKUPW(J,I,K)
!C3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
!C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
!C3A1----LOGARITHMIC MEAN INTERBLOCK TRANSMISSIVITY
               T2=HKUPW(J+1,I,K)
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
!C
!C3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
!C3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
               T2=HKUPW(J,I+1,K)
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
!C
!C4------RETURN
      RETURN
      END SUBROUTINE SGWF2UPW1HLOG

      SUBROUTINE SGWF2UPW1HUNCNF(K)
!C     ******************************************************************
!C-----COMPUTE CONSTANT PART OF HORIZONTAL CONDUCTANCE USING ARITHMETIC
!C-----MEAN CELL THICKNESS AND LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY.
!C-----THIS IS DIFFERENT FROM SGWF2UPW1HLOG FOR CONFINED LAYERS.
!C-----ACTIVATED BY LAYAVG=2
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC,BOTM,LBOTM
      USE GWFUPWMODULE,ONLY:HKUPW,CHANI,HANI,LAYTYPUPW
!C     ------------------------------------------------------------------
!C
      ZERO=0.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
      TWO=2.
!C
!C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
!C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
!C
!C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
!C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
!C
!C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         HYC1=HKUPW(J,I,K)
!C3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
!C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
!C3A1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
               HYC2=HKUPW(J+1,I,K)
               RATIO=HYC2/HYC1
               IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                  HYC=(HYC2-HYC1)/LOG(RATIO)
               ELSE
                  HYC=HALF*(HYC1+HYC2)
               END IF
               CR(J,I,K)=TWO*DELC(I)*HYC/(DELR(J+1)+DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
            CR(J,I,K)=ZERO
         END IF
!C
!C3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
!C3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
!C3B1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
               HYC2=HKUPW(J,I+1,K)
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
               CC(J,I,K)=TWO*DELR(J)*HYC/(DELC(I+1)+DELC(I))
            ELSE
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
      RETURN
      END SUBROUTINE SGWF2UPW1HUNCNF

      SUBROUTINE SGWF2UPW1HHARMCON(K)
!C     ******************************************************************
!C      COMPUTE THE HORIZONTAL CONDUCTANCE FOR CONFINED CELLS BASED ON THE
!C      HARMONIC AVEARGE K FOR ADJACENT CELLS.
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, iout, &
                       CC, CR, BOTM, LBOTM
      USE GWFUPWMODULE
!C     ------------------------------------------------------------------
!C
      ZERO=0.
      TWO=2.
!C1A-----Put cell thickness into CC.
       DO 200 I=1,NROW
       DO 200 J=1,NCOL
       TTOP=BOTM(J,I,LBOTM(K)-1)
       BBOT=BOTM(J,I,LBOTM(K))
       THCK=TTOP-BBOT
       CC(J,I,K)=THCK
  200  CONTINUE
!C
!C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
!C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
!C
!C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
!C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
!C
!C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         T1=HKUPW(J,I,K)*CC(J,I,K)
!C3A-----IF THIS IS NOT THE LAST COLUMN (RIGHTMOST), CALCULATE
!C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
               T2=HKUPW(J+1,I,K)*CC(J+1,I,K)
               CR(J,I,K)=TWO*T2*T1*DELC(I)/(T1*DELR(J+1)+T2*DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
!C3B-----IF THIS IS THE LAST COLUMN, SET BRANCH CONDUCTANCE=0.
            CR(J,I,K)=ZERO
         END IF
!C
!C3C-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
!C3C-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
               T2=HKUPW(J,I+1,K)*CC(J,I+1,K)
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
!C3D-----IF THIS IS THE LAST ROW, SET BRANCH CONDUCTANCE=0.
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
!C
!C4------RETURN
      RETURN
      END SUBROUTINE SGWF2UPW1HHARMCON

      SUBROUTINE SGWF2UPW1HLOGCON(K)
!C     ******************************************************************
!C-----COMPUTE THE HORIZONTAL CONDUCTANCE FOR CONFINED CELLS BASED ON THE
!C-----LOGARITHMIC  MEAN HYDRAULIC CONDUCTIVITY -- ACTIVATED BY LAYAVG=1
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Delr, Delc, CR, CC, &
                      BOTM, LBOTM
      USE GWFUPWMODULE
!C     ------------------------------------------------------------------
!C
      ZERO=0.
      TWO=2.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
!C1A-----Put cell thickness into CC.
       DO 200 I=1,NROW
       DO 200 J=1,NCOL
       TTOP=BOTM(J,I,LBOTM(K)-1)
       BBOT=BOTM(J,I,LBOTM(K))
       THCK=TTOP-BBOT
       CC(J,I,K)=THCK
  200  CONTINUE
!C
!C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
!C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
!C
!C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
!C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
!C
!C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         T1=HKUPW(J,I,K)*CC(J,I,K)
!C3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
!C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
!C3A1----LOGARITHMIC MEAN INTERBLOCK TRANSMISSIVITY
               T2=HKUPW(J+1,I,K)*CC(J+1,I,K)
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
!C
!C3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
!C3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
               T2=HKUPW(J,I+1,K)*CC(J,I+1,K)
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
!C
!C4------RETURN
      RETURN
      END SUBROUTINE SGWF2UPW1HLOGCON

      SUBROUTINE SGWF2UPW1HUNCNFCON(K)
!C     ******************************************************************
!C-----COMPUTE HORIZONTAL CONDUCTANCE FOR CONFINED CELLS USING ARITHMETIC
!C-----MEAN CELL THICKNESS AND LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY.
!C-----ACTIVATED BY LAYAVG=2
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,      ONLY:NCOL,NROW,IBOUND,CR,CC,DELR,DELC,BOTM,LBOTM
      USE GWFUPWMODULE,ONLY:HKUPW,CHANI,HANI,LAYTYPUPW
!C     ------------------------------------------------------------------
!C
      ZERO=0.
      HALF=0.5
      FRAC1=1.005
      FRAC2=0.995
!C1A-----Put cell thickness into CC.
       DO 200 I=1,NROW
       DO 200 J=1,NCOL
       TTOP=BOTM(J,I,LBOTM(K)-1)
       BBOT=BOTM(J,I,LBOTM(K))
       THCK=TTOP-BBOT
       CC(J,I,K)=THCK
  200  CONTINUE
!C
!C1------FOR EACH CELL CALCULATE BRANCH CONDUCTANCES FROM THAT CELL
!C1------TO THE ONE ON THE RIGHT AND THE ONE IN FRONT.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
!C
!C2------IF CELL IS DRY OR HK=0., SET CONDUCTANCE EQUAL TO 0 AND GO ON
!C2------TO NEXT CELL.
      IF(IBOUND(J,I,K).EQ.0 .OR. HKUPW(J,I,K).EQ.ZERO) THEN
         CR(J,I,K)=ZERO
         CC(J,I,K)=ZERO
      ELSE
!C
!C3------CELL IS WET -- CALCULATE TRANSMISSIVITY OF CELL.
         HYC1=HKUPW(J,I,K)
!C3A-----IF THIS IS NOT THE LAST COLUMN(RIGHTMOST) THEN CALCULATE
!C3A-----BRANCH CONDUCTANCE IN THE ROW DIRECTION (CR) TO THE RIGHT.
         IF(J.NE.NCOL) THEN
            IF(IBOUND(J+1,I,K).NE.0) THEN
!C3A1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
               HYC2=HKUPW(J+1,I,K)
               RATIO=HYC2/HYC1
               IF(RATIO.GT.FRAC1 .OR. RATIO.LT.FRAC2) THEN
                  HYC=(HYC2-HYC1)/LOG(RATIO)
               ELSE
                  HYC=HALF*(HYC1+HYC2)
               END IF
!C3A2----MULTIPLY LOGARITHMIC K BY ARITMETIC SATURATED THICKNESS.
               CR(J,I,K)=DELC(I)*HYC*(CC(J,I,K)+CC(J+1,I,K))/  &
                           (DELR(J+1)+DELR(J))
            ELSE
               CR(J,I,K)=ZERO
            END IF
         ELSE
            CR(J,I,K)=ZERO
         END IF
!C
!C3B-----IF THIS IS NOT THE LAST ROW (FRONTMOST) THEN CALCULATE
!C3B-----BRANCH CONDUCTANCE IN THE COLUMN DIRECTION (CC) TO THE FRONT.
         IF(I.NE.NROW) THEN
            IF(IBOUND(J,I+1,K).NE.0) THEN
!C3B1----LOGARITHMIC MEAN HYDRAULIC CONDUCTIVITY
               HYC2=HKUPW(J,I+1,K)
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
!C3B2----MULTIPLY LOGARITHMIC K BY ARITMETIC SATURATED THICKNESS.
               CC(J,I,K)=DELR(J)*HYC*(CC(J,I,K)+CC(J,I+1,K))/ &
                           (DELC(I+1)+DELC(I))
            ELSE
               CC(J,I,K)=ZERO
            END IF
         ELSE
            CC(J,I,K)=ZERO
         END IF
      END IF
  100 CONTINUE
!C
!C4------RETURN.
      RETURN
      END SUBROUTINE SGWF2UPW1HUNCNFCON

      SUBROUTINE SGWF2UPW1VCOND(K)
!C     ******************************************************************
!C     COMPUTE VERTICAL BRANCH CONDUCTANCE BETWEEN A LAYER AND THE NEXT
!C     LOWER LAYER FROM VERTICAL HYDRAULIC CONDUCTIVITY.
!C     ******************************************************************
!C
!C      SPECIFICATIONS:
!C     ------------------------------------------------------------------
      USE GLOBAL,        ONLY:NCOL,NROW,NLAY,IBOUND,HNEW,DELR,DELC, &
                              BOTM,LBOTM,LAYCBD,IOUT,STRT,CV
      USE GWFUPWMODULE
      use SimPHMFModule, only: ustop
!C
      DOUBLE PRECISION BBOT,TTOP,HHD
!C     ------------------------------------------------------------------
!C
      IF(K.EQ.NLAY) RETURN
      ZERO=0.
      HALF=0.5
!C
!C1------LOOP THROUGH ALL CELLS IN THE LAYER.
      DO 100 I=1,NROW
      DO 100 J=1,NCOL
      CV(J,I,K)=ZERO
      IF(IBOUND(J,I,K).NE.0 .AND. IBOUND(J,I,K+1).NE.0) THEN
!C
!C2------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL.
         IF(LAYVKAUPW(K).EQ.0) THEN
            HYC1=VKAUPW(J,I,K)
         ELSE
            HYC1=HKUPW(J,I,K)/VKAUPW(J,I,K)
         END IF
         IF(HYC1.GT.ZERO) THEN
!C3------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CELL BELOW.
            IF(LAYVKAUPW(K+1).EQ.0) THEN
               HYC2=VKAUPW(J,I,K+1)
            ELSE
               HYC2=(HKUPW(J,I,K+1)/VKAUPW(J,I,K+1))
            END IF
            IF(HYC2.GT.ZERO) THEN
!C
!C4------CALCULATE INVERSE LEAKANCE FOR CELL.  ICONCV FLAG PREVENTS
!C4------CV FROM BEING HEAD DEPENDENT.
               BBOT=BOTM(J,I,LBOTM(K))
               TTOP=BOTM(J,I,LBOTM(K)-1)
               IF(LAYSTRT(K).NE.0) TTOP=STRT(J,I,K)
               IF(LAYTYPUPW(K).NE.0 .AND. ICONCV.EQ.0) THEN
                  HHD=HNEW(J,I,K)
!                  IF(HHD.LT.TTOP) TTOP=HHD   !RGN 6/23/09
               END IF
               BOVK1=(TTOP-BBOT)*HALF/HYC1
!C
!C5------CALCULATE INVERSE LEAKANCE FOR CELL BELOW.
               BBOT=BOTM(J,I,LBOTM(K+1))
               TTOP=BOTM(J,I,LBOTM(K+1)-1)
               IF(LAYSTRT(K+1).NE.0) TTOP=STRT(J,I,K+1)
               BBB=(TTOP-BBOT)*HALF
!C
!C5A-----IF CELL BELOW IS NOT SATURATED, DO NOT INCLUDE ITS CONDUCTANCE
!C5A-----IN THE VERTICAL CONDUCTANCE CALULATION, EXCEPT THAT THE NOCVCO
!C5A-----AND ICONCV FLAGS TURN OFF THIS CORRECTION.
               IF(LAYTYPUPW(K+1).NE.0 &
                         .AND.NOCVCO.EQ.0 .AND. ICONCV.EQ.0) THEN
                  HHD=HNEW(J,I,K+1)
!                  IF(HHD.LT.TTOP) BBB=ZERO   !RGN 6/23/09
               END IF
               BOVK2=BBB/HYC2
!C
!C6------CALCULATE VERTICAL HYDRAULIC CONDUCTIVITY FOR CONFINING BED.
               IF(LAYCBD(K).NE.0) THEN
                  IF(VKCB(J,I,LAYCBD(K)).GT.ZERO) THEN
!C
!C7------CALCULATE INVERSE LEAKANCE FOR CONFINING BED.
                     BBB=BOTM(J,I,LBOTM(K))-BOTM(J,I,LBOTM(K)+1)
                     IF(BBB.LT.ZERO) THEN
                        WRITE(IOUT,45) K,I,J
   45                   FORMAT(1X,/1X, &
       'Negative confining bed thickness below cell (Layer,row,col)', &
                       I4,',',I5,',',I5)
            WRITE(IOUT,46) BOTM(J,I,LBOTM(K)),BOTM(J,I,LBOTM(K)+1)
   46       FORMAT(1X,'Top elevation, bottom elevation:',1P,2G13.5)
                        CALL USTOP(' ')
                     END IF
                     CBBOVK=BBB/VKCB(J,I,LAYCBD(K))
                     CV(J,I,K)=DELR(J)*DELC(I)/(BOVK1+CBBOVK+BOVK2)
                  END IF
               ELSE
                  CV(J,I,K)=DELR(J)*DELC(I)/(BOVK1+BOVK2)
               END IF
            END IF
         END IF
      END IF
  100 CONTINUE
!C
!C8------RETURN.
      RETURN
      END SUBROUTINE SGWF2UPW1VCOND

      SUBROUTINE SGWF2UPW1PNT(Igrid)
      USE GWFUPWMODULE
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER Igrid
!     ------------------------------------------------------------------
! Cell property data
        Sn=>Gwfupwdat(IGRID)%Sn
        So=>Gwfupwdat(IGRID)%So
        IUPWCB=>Gwfupwdat(IGRID)%IUPWCB
        IWDFLG=>Gwfupwdat(IGRID)%IWDFLG
        IWETIT=>Gwfupwdat(IGRID)%IWETIT
        IHDWET=>Gwfupwdat(IGRID)%IHDWET
        IPHDRY=>Gwfupwdat(IGRID)%IPHDRY
        ISFAC=>Gwfupwdat(IGRID)%ISFAC
        ICONCV=>Gwfupwdat(IGRID)%ICONCV
        ITHFLG=>Gwfupwdat(IGRID)%ITHFLG
        NOCVCO=>Gwfupwdat(IGRID)%NOCVCO
        NOVFC=>Gwfupwdat(IGRID)%NOVFC
        WETFCT=>Gwfupwdat(IGRID)%WETFCT
        LAYTYPUPW=>Gwfupwdat(IGRID)%LAYTYPUPW
        LAYAVG=>Gwfupwdat(IGRID)%LAYAVG
        CHANI=>Gwfupwdat(IGRID)%CHANI
        LAYVKAUPW=>Gwfupwdat(IGRID)%LAYVKAUPW
        LAYWET=>Gwfupwdat(IGRID)%LAYWET
        LAYSTRT=>Gwfupwdat(IGRID)%LAYSTRT
        LAYFLG=>Gwfupwdat(IGRID)%LAYFLG
        VKAUPW=>Gwfupwdat(IGRID)%VKAUPW
        VKCB=>Gwfupwdat(IGRID)%VKCB
        SC1=>Gwfupwdat(IGRID)%SC1
        SC2UPW=>Gwfupwdat(IGRID)%SC2UPW
        HANI=>Gwfupwdat(IGRID)%HANI
        WETDRY=>Gwfupwdat(IGRID)%WETDRY
        HKUPW=>Gwfupwdat(IGRID)%HKUPW
        IBOUND2=>Gwfupwdat(IGRID)%IBOUND2
      END SUBROUTINE SGWF2UPW1PNT

      SUBROUTINE SGWF2UPW1PSV(Igrid)
      USE GWFUPWMODULE
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER Igrid
!     ------------------------------------------------------------------
! Cell property data
        Gwfupwdat(IGRID)%Sn=>Sn
        Gwfupwdat(IGRID)%So=>So
        Gwfupwdat(IGRID)%IUPWCB=>IUPWCB
        Gwfupwdat(IGRID)%IWDFLG=>IWDFLG
        Gwfupwdat(IGRID)%IWETIT=>IWETIT
        Gwfupwdat(IGRID)%IHDWET=>IHDWET
        Gwfupwdat(IGRID)%IPHDRY=>IPHDRY
        Gwfupwdat(IGRID)%ISFAC=>ISFAC
        Gwfupwdat(IGRID)%ICONCV=>ICONCV
        Gwfupwdat(IGRID)%ITHFLG=>ITHFLG
        Gwfupwdat(IGRID)%NOCVCO=>NOCVCO
        Gwfupwdat(IGRID)%NOVFC=>NOVFC
        Gwfupwdat(IGRID)%WETFCT=>WETFCT
        Gwfupwdat(IGRID)%LAYTYPUPW=>LAYTYPUPW
        Gwfupwdat(IGRID)%LAYAVG=>LAYAVG
        Gwfupwdat(IGRID)%CHANI=>CHANI
        Gwfupwdat(IGRID)%LAYVKAUPW=>LAYVKAUPW
        Gwfupwdat(IGRID)%LAYWET=>LAYWET
        Gwfupwdat(IGRID)%LAYSTRT=>LAYSTRT
        Gwfupwdat(IGRID)%LAYFLG=>LAYFLG
        Gwfupwdat(IGRID)%VKAUPW=>VKAUPW
        Gwfupwdat(IGRID)%VKCB=>VKCB
        Gwfupwdat(IGRID)%SC1=>SC1
        Gwfupwdat(IGRID)%SC2UPW=>SC2UPW
        Gwfupwdat(IGRID)%HANI=>HANI
        Gwfupwdat(IGRID)%WETDRY=>WETDRY
        Gwfupwdat(IGRID)%HKUPW=>HKUPW
        Gwfupwdat(IGRID)%IBOUND2=>IBOUND2
!
      END SUBROUTINE SGWF2UPW1PSV

end module UpwSubsModule

module unused
!      SUBROUTINE GWF2UPWFMS(KITER,KSTP,KPER,IGRID)
!C     ******************************************************************
!C     ADD LEAKAGE CORRECTION AND STORAGE TO HCOF AND RHS, AND CALCULATE
!C     CONDUCTANCE AS REQUIRED.
!C     ******************************************************************
!C
!C        SPECIFICATIONS:
!C     ------------------------------------------------------------------
!      USE GLOBAL,      ONLY:NCOL,NROW,NLAY,IBOUND,BOTM,NBOTM,DELR,DELC,
!     1                      LBOTM,CV,HNEW,RHS,HCOF,HOLD,ISSFLG,IOUT
!      USE GWFBASMODULE,ONLY:DELT
!      USE GWFNWTMODULE,ONLY:A, ICELL, IA
!      USE GWFUPWMODULE
!      DOUBLE PRECISION, EXTERNAL :: DHORIZUPW
!      DOUBLE PRECISION ZERO
!      DOUBLE PRECISION HTMP, TP, BT, TLED, ONE, SOLD, SNEW, STRG
!      DOUBLE PRECISION RHO1, RHO2, HLD, THICK, dS
!C     ------------------------------------------------------------------
!C
!C1------SET POINTERS TO DATA, GET STEADY-STATE FLAG FOR STRESS PERIOD,
!C1------DEFINE CONSTANT.
!      CALL SGWF2UPW1PNT(IGRID)
!      ISS=ISSFLG(KPER)
!      ONE=1.0D0
!      ZERO = 0.0D0
!C
!C2------IF THE STRESS PERIOD IS TRANSIENT, ADD STORAGE TO HCOF AND RHS
!      IF(ISS.EQ.0) THEN
!         TLED=ONE/DBLE(DELT)
!         DO 200 K=1,NLAY
!C3------CHECK OLD AND NEW HEADS TO DETERMINE
!C3------WHEN TO USE PRIMARY AND SECONDARY STORAGE
!            DO 180 I=1,NROW
!            DO 180 J=1,NCOL
!C
!C4A-----IF THE CELL IS EXTERNAL THEN SKIP IT.
!            IF(IBOUND(J,I,K).LE.0) GO TO 180
!            TP=dble(BOTM(J,I,LBOTM(K)-1))
!            BT=dble(BOTM(J,I,LBOTM(K)))
!            THICK = (TP-BT)
!            HTMP = HNEW(J,I,K)
!            HLD = DBLE(HOLD(J,I,K))
!            RHO2 = SC2UPW(J,I,K)*TLED
!!            IF ( HTMP.GT.TP ) RHO2 = SC1(J,I,K)*TLED
!            RHO1 = SC1(J,I,K)*TLED
!!            IF ( HLD.GT.TP ) RHO1 = SC1(J,I,K)*TLED
!            dS = DHORIZUPW(HTMP, TP, BT, K)
!
!! Add derivatives to jacobian
!            ij = Icell(J,I,K)
!! Derivative for HCOF
!            A(IA(ij)) = A(IA(ij)) - dS*THICK*RHO1*HTMP
!! Derivative for RHS
!            A(IA(ij)) = A(IA(ij)) - THICK*RHO2*dS
!     +                            + THICK*RHO1*dS*HLD
!
!C
!C6------ADD STORAGE TERMS TO RHS AND HCOF.
!        HCOF(J,I,K) = HCOF(J,I,K) - Sn(ij)*THICK*RHO1
!        RHS(J,I,K) = RHS(J,I,K) + THICK*RHO2*(Sn(ij)-So(ij)) -
!     +                            Sn(ij)*THICK*RHO1*HLD
!C
!  180       CONTINUE
!C
!  200    CONTINUE
!      END IF
!C
!C10-----RETURN
!      RETURN
!      END
!C
!C
!!     -----------------------------------------------------------------
!!     Updates saturation for previous time step.
!
!      SUBROUTINE GWF2UPW1AD(IGRID)
!      USE GWFUPWMODULE, ONLY: Sn, So
!      USE GWFNWTMODULE, ONLY: Numactive
!      IMPLICIT NONE
!!     -----------------------------------------------------------------
!!     ARGUMENTS
!      INTEGER IGRID
!!     ------------------------------------------------------------------
!
!!     LOCAL VARIABLES
!!     -----------------------------------------------------------------
!      INTEGER ij
!!     -----------------------------------------------------------------
!C
!C------SET POINTERS FOR THE CURRENT GRID.
!      CALL SGWF2UPW1PNT(IGRID)
!C------Set old saturation to new saturation.
!      DO ij = 1, Numactive
!        So(ij) = Sn(ij)
!      END DO
!      RETURN
!      END SUBROUTINE
!!
!!
!!     SUBROUTINE GWF2UPWUPDATE. UPDATE VALUES AFTER OUTER ITERATION.
!      SUBROUTINE GWF2UPWUPDATE(Itest, Igrid)
!      USE GLOBAL, ONLY: Ncol, Nrow, Nlay, Ibound, HNEW
!      USE GWFNWTMODULE, ONLY: A, IA, Numactive, Diag, HITER
!      IMPLICIT NONE
!!     ------------------------------------------------------------------
!!     SPECIFICATIONS:
!!     ------------------------------------------------------------------
!      EXTERNAL SGWF2UPW1PNT
!      DOUBLE PRECISION, EXTERNAL:: Sat_thick
!!     ------------------------------------------------------------------
!!     ARGUMENTS
!!     ------------------------------------------------------------------
!      INTEGER Igrid, ij, ic, ir, il, Itest, i, i1, i2
!!     -----------------------------------------------------------------
!C------SET POINTERS FOR THE CURRENT GRID.
!      CALL SGWF2UPW1PNT(Igrid)
!      DO ij = 1, Numactive
!        il = Diag(ij, 1)
!        ir = Diag(ij, 2)
!        ic = Diag(ij, 3)
!        IF (Itest==1) HNEW(ic,ir,il) = HITER(ic,ir,il)
!        I1 = IA(ij)
!        I2 = IA(ij+1)-1
!        DO i = I1, I2
!          A(i) = 0.0D0
!        END DO
!      END DO
!      CALL Sn_update()
!      END SUBROUTINE
!!     -----------------------------------------------------------------
!!     Updates saturation for latest iteration.
!
!      SUBROUTINE Sn_update()
!      USE GWFUPWMODULE
!      USE GWFNWTMODULE, ONLY: Diag, Numactive
!      USE GLOBAL,      ONLY: Iout,HNEW,BOTM,LBOTM,NLAY
!      IMPLICIT NONE
!!     ------------------------------------------------------------------
!!     SPECIFICATIONS:
!!     -----------------------------------------------------------------
!      DOUBLE PRECISION, EXTERNAL :: Sat_thick
!!     -----------------------------------------------------------------
!!     LOCAL VARIABLES
!!     -----------------------------------------------------------------
!      DOUBLE PRECISION HH, TP, BT
!      INTEGER ij, ic, ir, il
!!     -----------------------------------------------------------------
!      DO ij = 1, Numactive
!        il = Diag(ij, 1)
!        ir = Diag(ij, 2)
!        ic = Diag(ij, 3)
!        HH = HNEW(ic,ir,il)
!        TP = BOTM(ic,ir,LBOTM(il)-1)
!        BT = BOTM(ic,ir,LBOTM(il))
!        Sn(ij) = Sat_thick(HH, TP, BT, il)
!      END DO
!      RETURN
!      END SUBROUTINE
!!
!!     ------------------------------------------------------------------
!!
!!
!!     -----------------------------------------------------------------
!!
!      FUNCTION DHORIZUPW(Hup, Ttop, Bbot, il)
!! RETURNS DERIVATIVE OF HORIZONTAL CONDUCTANCE BASED ON SMOOTH FUNCTION
!! FUNCTION IS CALCULATED IN UPW PACKAGE IN SUBROUTINE SAT_THICK
!      USE GWFNWTMODULE, ONLY: Thickfact
!      USE GWFUPWMODULE, ONLY: LAYTYPUPW
!      IMPLICIT NONE
!!     ------------------------------------------------------------------
!!     SPECIFICATIONS:
!!     ------------------------------------------------------------------
!!     ------------------------------------------------------------------
!!     ARGUMENTS
!!     -----------------------------------------------------------------
!      DOUBLE PRECISION Hup, Ttop, Bbot
!!     -----------------------------------------------------------------
!!     LOCAL VARIABLES
!!     -----------------------------------------------------------------
!      DOUBLE PRECISION factor, x, s, v, cof1, cof2, EPS, ACOF, Y
!      DOUBLE PRECISION EPSQD, z
!      DOUBLE PRECISION DHORIZUPW
!      INTEGER il
!!     -----------------------------------------------------------------
!C-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
!      DHORIZUPW = 0.0D0
!      IF ( LAYTYPUPW(il).LE.0 ) RETURN
!      EPS = Thickfact
!      ACOF = 1.0 / (1.0 - EPS)
!      x = (Hup-bbot)/(TTOP-BBOT)
!      IF ( x.LT.1.0d-9 )  x = 1.0d-9
!      IF(X.LT.EPS)THEN
!        Y = ACOF * X / (EPS*(Ttop-Bbot))
!      ELSEIF(X.LT.1.0D0-EPS)THEN
!        Y = ACOF /(Ttop-Bbot)
!      ELSEIF(X.LT.1.0D0)THEN
!        X = 1.0 - X
!        Y = - ACOF * x / (EPS * (Ttop - Bbot))
!        Y = 1.0-Y
!      ELSE
!        Y = 0.0
!      ENDIF
!      factor = Y
!      DHORIZUPW = factor
!      END FUNCTION DHORIZUPW
!!
!!
!!     ------------------------------------------------------------------
!!
!      SUBROUTINE GWF2UPW1DA(Igrid)
!      USE GWFUPWMODULE
!      IMPLICIT NONE
!!     ------------------------------------------------------------------
!!     SPECIFICATIONS:
!!     ------------------------------------------------------------------
!!     ARGUMENTS
!!     ------------------------------------------------------------------
!      INTEGER Igrid
!!     ------------------------------------------------------------------
!! Deallocate UPW data.
!        DEALLOCATE(Gwfupwdat(IGRID)%Sn)
!        DEALLOCATE(Gwfupwdat(IGRID)%So)
!        DEALLOCATE(Gwfupwdat(IGRID)%IUPWCB)
!        DEALLOCATE(Gwfupwdat(IGRID)%IWDFLG)
!        DEALLOCATE(Gwfupwdat(IGRID)%IWETIT)
!        DEALLOCATE(Gwfupwdat(IGRID)%IHDWET)
!        DEALLOCATE(Gwfupwdat(IGRID)%IPHDRY)
!        DEALLOCATE(Gwfupwdat(IGRID)%ISFAC)
!        DEALLOCATE(Gwfupwdat(IGRID)%ICONCV)
!        DEALLOCATE(Gwfupwdat(IGRID)%ITHFLG)
!        DEALLOCATE(Gwfupwdat(IGRID)%NOCVCO)
!        DEALLOCATE(Gwfupwdat(IGRID)%NOVFC)
!        DEALLOCATE(Gwfupwdat(IGRID)%WETFCT)
!        DEALLOCATE(Gwfupwdat(IGRID)%LAYTYPUPW)
!        DEALLOCATE(Gwfupwdat(IGRID)%LAYAVG)
!        DEALLOCATE(Gwfupwdat(IGRID)%CHANI)
!        DEALLOCATE(Gwfupwdat(IGRID)%LAYVKAUPW)
!        DEALLOCATE(Gwfupwdat(IGRID)%LAYWET)
!        DEALLOCATE(Gwfupwdat(IGRID)%LAYSTRT)
!        DEALLOCATE(Gwfupwdat(IGRID)%LAYFLG)
!        DEALLOCATE(Gwfupwdat(IGRID)%VKAUPW)
!        DEALLOCATE(Gwfupwdat(IGRID)%VKCB)
!        DEALLOCATE(Gwfupwdat(IGRID)%SC1)
!        DEALLOCATE(Gwfupwdat(IGRID)%SC2UPW)
!        DEALLOCATE(Gwfupwdat(IGRID)%HANI)
!        DEALLOCATE(Gwfupwdat(IGRID)%WETDRY)
!        DEALLOCATE(Gwfupwdat(IGRID)%HKUPW)
!        DEALLOCATE(Gwfupwdat(IGRID)%IBOUND2)
!      END SUBROUTINE GWF2UPW1DA

!C4------RETURN.
!
!      DOUBLE PRECISION FUNCTION SAT_THICK(Hup,Ttop,Bbot,il)
!! RETURNS SATURATED THICKNESS OF CELL BASED ON SMOOTH FUNCTION
!      USE GWFUPWMODULE
!      USE GWFNWTMODULE, ONLY: Thickfact
!      USE GLOBAL, ONLY: IOUT
!      IMPLICIT NONE
!!     ------------------------------------------------------------------
!!     SPECIFICATIONS:
!!     ------------------------------------------------------------------
!!     LOCAL VARIABLES
!!     -----------------------------------------------------------------
!      INTEGER ic, ir, il, iltyp, METHOD1 ! make METHOD global
!      DOUBLE PRECISION hup, bbot, zero, ttop, factor, x, EPS, ACOF
!      DOUBLE PRECISION cof1, cof2, s, v, factor1, factor2, Y, Z, EPSQD
!!     -----------------------------------------------------------------
!! Calculate saturated thickenss
!      zero = 0.0D0
!!      v = (ttop-bbot)
!      SAT_THICK = 1.0D0
!      IF ( LAYTYPUPW(il).LE.0 ) RETURN
!C-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
!      EPS = Thickfact
!      ACOF = 1.0 / (1.0 - EPS)
!      x = (Hup-bbot)/(TTOP-BBOT)
!      IF ( x.LT.1.0d-9 )  x = 1.0d-9
!      IF(X.LT.EPS)THEN
!        Y = ACOF *0.5/EPS * X**2
!      ELSEIF(X.LT.1.0-EPS)THEN
!        Y = ACOF * X + (1.0-ACOF)*0.5
!      ELSEIF(X.LT.1.0)THEN
!        X = 1.0 - X
!        Y = ACOF *0.5/EPS * X**2
!        Y = 1.0-Y
!      ELSE
!        Y = 1.0
!      ENDIF
!      factor = Y
!!
!      SAT_THICK = factor
!      END FUNCTION SAT_THICK
!!
!!     -----------------------------------------------------------------
!
    end module unused
