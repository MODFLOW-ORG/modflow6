 
!
!-------SUBROUTINE GWF2NWT1AR
!
      SUBROUTINE GWF2NWT1AR(In, Mxiter, Iunitlak, Igrid, ibotavg)
!
!------NEWTON SOLVER VERSION NUMBER 1.0.9:  July 1, 2014
!      RICHARD G. NISWONGER
      USE GLOBAL,     ONLY:NCOL,NROW,NLAY,ITRSS,LAYHDT,LAYHDS,LAYCBD,
     1                     NCNFBD,IBOUND,BUFF,BOTM,NBOTM,DELR,DELC,IOUT,
     2                     LBOTM,HNEW
      USE GWFBASMODULE,ONLY:HDRY
      USE GWFNWTMODULE
      use utl7module, only: urdcom, urword
      use SimPHMFModule, only: ustop
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     ------------------------------------------------------------------
      INTRINSIC INT
!      EXTERNAL URDCOM, URWORD
!      EXTERNAL SGWF2NWT1PSV
!     ------------------------------------------------------------------
!     ARGUMENTS
!     ------------------------------------------------------------------
      INTEGER In, Igrid, Mxiter, L, NRC, Iunitlak
      integer, intent(out) :: ibotavg
!     ------------------------------------------------------------------
!     LOCAL VARIABLES
!     ------------------------------------------------------------------
      INTEGER lloc, LLOCSAVE, istart, istop, i, ic, ir, il, jj
      CHARACTER(LEN=300) line
      REAL relaxdum
      REAL ZERO
      INTEGER IANAME,KHANI,N,KK,nc,nr,nl,j,k,NCNVRT,NHANI,NWETD
! Memory use variables
      INTEGER lrwrk,liwrk,NODES,MBLACK,NJAF
      REAL Memuse1,Memuse2
      double precision :: akappadum, r, toldum, ftoldum, Thickdum, 
     &                    gammadum, amomentdum, thetadum, Breducdum, 
     &                    Btoldum
!     ------------------------------------------------------------------
!

      CALL URDCOM(In, Iout, line)
      lloc = 1
      ALLOCATE (Tol, Ftol, RMS2, RMS1, Iierr,IFDPARAM,ICNVGFLG)
      ALLOCATE (ITER1,THETA,THICKFACT,BTOL,Numtrack)
      ALLOCATE (RMSAVE)
      ALLOCATE (Numnonzero, Numactive, Numcell, II)
      ALLOCATE (Akappa,Gamma,Amomentum,Btrack,Breduc,Numtrack)
      ALLOCATE (Nonmeth, Linmeth, IPRNWT, Itreal, Ibt)
      ALLOCATE (IBOTAV)
      allocate (MxIterInner, HCloseLinear, RCloseLinear)
      allocate (iaclNwt, norderNwt, levelNwt, northNwt, MaxBackIterNwt, 
     &          IRedSysNwt, IDropTolNwt, BackFlagNwt, BackTolNwt,
     &          BackReduceNwt, EpsrnNwt, RrctolNwt)
      ! for MF5to6, assign defaults to new variables
      MxIterInner = 50
      HCloseLinear = 1.0d-4
      RCloseLinear = 1.0d-10
!1------IDENTIFY PACKAGE AND INITIALIZE.
      write(*,*)'Processing NWT package input...'
      WRITE (Iout, 9001) In
 9001 FORMAT (1X, /' NWT1 -- Newton Solver, ',
     +    'VERSION 1.0.9, 07/01/2014', /, 9X, 'INPUT READ FROM UNIT',
     +        I3,/)
      i = 1
      Itreal = 0
      Ibt = 0
      RMS2 = 0.0D0
      RMS1 = 0.0D0
      RMSAVE = 0.0D0
      Numactive = 0
      Numcell = 0
      akappadum = 0.0   
      toldum = 1.0e-4
      ftoldum = 100.0
      Mxiter = 100
      Thickdum = 1.0e-4
      Linmeth = 2     ! Linmeth=1 GMRES; Linmeth=2 XMD; Linmeth=3 SAMG 
      IPRNWT = 1     ! Iteration stats (>0 prints)
      IBOTAV = 1     ! Doesn't reset to bottom
      LLOC=1
      Numtrack = 0
      Btoldum = 1.0
      Breducdum = 1.0
      iaclNwt = 2
      norderNwt = 1
      levelNwt = 0
      northNwt = 2
      MaxBackIterNwt = 1 
      IRedSysNwt = 0
      IDropTolNwt = 1
      BackFlagNwt = 0
      BackTolNwt = 1.0d0
      BackReduceNwt = 0.05d0
      EpsrnNwt = 1.0d-3
      RrctolNwt = 0.0d0
      ICNVGFLG = 0
      CALL URWORD(line, lloc, istart, istop, 3, i, toldum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, ftoldum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, Mxiter, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, Thickdum, Iout, In)
!      CALL URWORD(line, lloc, istart, istop, 2, Nonmeth, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, Linmeth, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, IPRNWT, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, IBOTAV, r, Iout, In)
      ibotavg = IBOTAV
      ITER1 = Mxiter
C
C3B-----GET OPTIONS.
      IFDPARAM=0
   20 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'SIMPLE') THEN
         IFDPARAM=1
         WRITE(IOUT,21)
   21    FORMAT(1X,'SIMPLE OPTION:',/,
     1     1X,'DEFAULT SOLVER INPUT VALUES REFLECT NEARLY LINEAR MODEL')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'MODERATE') THEN
         IFDPARAM=2
         WRITE(IOUT,23)
   23    FORMAT(1X,'MODERATE OPTION:',/,1X,'DEFAULT SOLVER',
     1         ' INPUT VALUES REFLECT MODERETELY NONLINEAR MODEL')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'COMPLEX') THEN
         IFDPARAM=3
         WRITE(IOUT,25)
   25    FORMAT(1X,'COMPLEX OPTION:',/,1X,'DEFAULT SOLVER',
     1 ' INPUT VALUES REFLECT STRONGLY NONLINEAR MODEL')
      ELSE IF(LINE(ISTART:ISTOP).EQ.'SPECIFIED') THEN
         IFDPARAM=4
         WRITE(IOUT,26)
   26    FORMAT(1X,'SPECIFIED OPTION:',/,1X,'SOLVER INPUT',
     1 ' VALUES ARE SPECIFIED BY USER')
      END IF
      LLOCSAVE = LLOC
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'CONTINUE') THEN
        ICNVGFLG = 1
      ELSE
        LLOC = LLOCSAVE
      END IF
        
      
!      IF(LLOC.LT.200) GO TO 20    
!
! Don't need to read these when using default Options.
      IF ( IFDPARAM.EQ.4 ) THEN
      CALL URWORD(line, lloc, istart, istop, 3, i, thetadum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, akappadum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, gammadum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, amomentdum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 2, Btrack, r, Iout, In)
      IF ( BTRACK.GT.0 ) THEN
      CALL URWORD(line, lloc, istart, istop, 2, Numtrack, r, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, Btoldum, Iout, In)
      CALL URWORD(line, lloc, istart, istop, 3, i, Breducdum, Iout, In)
      END IF
      ELSEIF ( IFDPARAM.EQ.1 ) THEN
! Set values based on default and Option keyword.
        thetadum = 0.97
        akappadum = 0.0001
        gammadum = 0.0
        amomentdum = 0.0
        Btrack = 0
        Numtrack = 20
        Btoldum = 1.5
        Breducdum = 0.97
      ELSEIF ( IFDPARAM.EQ.2 ) THEN
        thetadum = 0.90
        akappadum = 0.0001
        gammadum = 0.00
        amomentdum = 0.1
        Btrack = 0
        Numtrack = 20
        Btoldum = 1.1
        Breducdum = 0.9
      ELSEIF ( IFDPARAM.EQ.3 ) THEN
        thetadum = 0.85
        akappadum = 0.00001
        gammadum = 0.0
        amomentdum = 0.1
        Btrack = 1
        Numtrack = 50
        Btoldum = 1.1
        Breducdum = 0.7
      ELSE
        Write(iout,*)
        Write(iout,*)'***Erroneous value for Input value "Options."***'
        Write(iout,*)'Check input. Model Stopping.'
        Write(iout,*) 
        CALL USTOP(' ')
      END IF
      !
      ! Assign backtracking-related variables to be used by SMS
      if (Btrack > 0) then
        BackFlagNwt = Btrack
        MaxBackIterNwt = Numtrack
        BackTolNwt = Btoldum
        BackReduceNwt = Breducdum
      endif
 !     
 !     IF ( Nonmeth==1 )Then
 !       Write(iout,*) '***Newton Linearization will be used***'
 !       Write(iout,*)
 !     ELSEIF ( Nonmeth==0 )Then
 !       Write(iout,*) '***Picard Linearization will be used***'
 !       Write(iout,*)
 !     ELSE
 !       Write(iout,*) '***Incorrect value for variable Nonmeth was ',
 !    +                'specified. Check input.***'
 !       Write(iout,*)
 !       Call USTOP('  ')
 !     END IF
      Nonmeth = 1
!
      IF ( Linmeth==1 )Then
        Write(iout,*) '***GMRES linear solver will be used***'
        Write(iout,*)
      ELSEIF ( Linmeth==2 )Then
        Write(iout,*) '***XMD linear solver will be used***'
        Write(iout,*)
       ELSEIF ( Linmeth==3 )Then
        Write(iout,*) '***SAMG linear solver will be used***'
        Write(iout,*)
      ELSE
        Write(iout,*) '***Incorrect value for Linear solution method ',
     +                'specified. Check input.***'
        Write(iout,*)
        Call USTOP('  ')
      END IF
!
      Thickfact = Thickdum
      Btol = Btoldum
      Breduc = Breducdum
      Theta = Thetadum
      Akappa = akappadum
      Gamma = gammadum
      Amomentum = amomentdum
      IF ( Theta.LT.CLOSEZERO ) Theta = 0.9D0
      Tol = toldum
      Ftol = ftoldum
!2--ECHO NWT INPUT
      WRITE(IOUT,9010) Tol,Ftol,MXITER
      WRITE(IOUT,9011) THETA,AKAPPA,GAMMADUM,AMOMENTUM 
 9010 FORMAT(1X,'  CONVERGENCE CRITERION OF',E15.6,' FOR HEAD SOLUTION',
     +      /1X,'  AND A TOLERANCE OF',E15.6,' FOR FLOW SOLUTION AND ',
     +      /1X,'  A MAXIMUM OF ',I5,' OUTER ITERATIONS. ',//)
 9011 FORMAT(1X,'  D-B-D REDUCTION FACTOR OF ',E15.6,' AND ',
     +      /1X,'  A D-B-D INCREASE FACTOR OF ',E15.6,' AND ',
     +      /1X,'  A D-B-D RELAXATION OF ',E15.6,' AND ', 
     +      /1X,'  A MOMENTUM FACTOR OF ',E15.6,' .',//)
      IF ( BTRACK.GT.0 ) THEN 
        WRITE(IOUT,9012) Numtrack,BTOL,BREDUC
      ELSE
        WRITE(IOUT,*)'***BACKTRACKING IS INACTIVE***'
      END IF
 9012 FORMAT(1X,'  BACKTRACKING IS ACTIVE ',
     +      /1X,'  THE MAXIMUM NUMBER OF BACKTRACKS IS ',I5,' AND ',
     +      /1X,'  THE BACKTRACKING TOLERANCE IS ',E15.6, ' AND',
     +      /1X,'  THE BACKTRACKING REDUCTION FACTOR IS ',E15.6,/)
!
!3-----ALLOCATE SPACE
      
      ALLOCATE (Icell(Ncol, Nrow, Nlay))
      ALLOCATE (Diag(Ncol*Nrow*Nlay, 3), II, Hiter(Ncol, Nrow, Nlay))
      Numnonzero = 0
      II = 0
      ICELL = 0
      DIAG = 0
 ! Check heads and set to be above bottom. 
      Do il = 1, Nlay
        Do ir = 1, Nrow
          Do ic = 1, Ncol
            IF ( IBOUND(ic,ir,il).GT.0 ) THEN
              IF ( dble(BOTM(ic,ir,LBOTM(il)-1)) - 
     +             dble(BOTM(ic,ir,LBOTM(il))).LT.100.0*Thickfact ) THEN
                WRITE(IOUT,*) 'Extremely thin cell for Column = ',ic,
     +                         ' and Row = ',ir,' and Layer = ',il,
     +                         ' Check input, Setting IBOUND = 0'
                IBOUND(ic,ir,il) = 0
              END IF    
 ! these next three lines could cause slow convergence when using a solution for IC.            
!              IF ( HNEW(ic,ir,il).LT.BOTM(ic,ir,LBOTM(il)) ) THEN
!                HNEW(ic,ir,il) = BOTM(ic,ir,LBOTM(il))+HEPS
!              END IF
            END IF
            Hiter(ic,ir,il) = HNEW(ic,ir,il)
          End do
        End do
      End do
 !  Determine the number of active cells and then numnber of elements
 !  in linear matrix for allocating arrays.
      CALL ORDERCELL()
      CALL COUNTACTIVE(jj)
      IF ( Numactive.LT.2 ) THEN
        WRITE(Iout,*)'MODFLOW-NWT does run with single-cell models. ',
     +               'Model Stopping.'
        CALL USTOP('')
      END IF
      Numnonzero = jj
      Numcell = Numactive
! Allocate global linear solver arrrays
! These are allocated exactly based on active cells
! and constant head cells for now. 
      NODES = NUMACTIVE
cmi
c     MBLACK = NODES
c     NJAF = 7 * Numnonzero
c     liwrk = 3*NODES + 3*mblack + njaf + 1    ! ldcomb = .false.
c     lrwrk = 7/2*mblack + 2*(14+1)*mblack + 14
cmi


      mblack = nodes
c     if (iredsys.eq.1) mblack = nodes * 0.5 + 1
c     njaf = 7 * nja
      njaf = 7 * numnonzero
c     liwrk = 3*nodes + 3*mblack + njaf + 1    ! ldcomb = .false.
c     if (ldcomb) liwrk = 3*nodes + 4*mblack + 2*njaf + 1  !        = .true.
      liwrk = 3*nodes + 4*mblack + 2*njaf + 1  !        = .true.

c     lrwrk = 3*mblack + 2*(north+1)*mblack + north
      lrwrk = 3*mblack + 2*(14+1)*mblack + 14

cmi
C
 !     memuse1 = 0
 !     memuse2 = 0
 !     memuse1 = (11*ncol*nrow*nlay+8*
 !    +           (Numnonzero+6*Numcell+6*Numcell)+
 !    +            4*((Numcell+1)+Numnonzero))
!      memuse1 = memuse1/1.0e9
 !     IF ( Linmeth==2 ) THEN
 !       memuse2 = 4
 !       memuse2 = 16*memuse2*NODES*7/1.0e9
 !     ELSE
 !       memuse2 = 8*(NJAF)/1.0e9
 !       memuse2 = memuse2 + 8*(LRWRK+NODES)/1.0e9
 !       memuse2 = memuse2 + 4*(NODES+LIWRK+30)/1.0e9
 !     END IF
!      Write(*,*) 'Gigabytes required for Nonlinear Solver= ',
!     +            memuse1
!      Write(*,*) 'Gigabytes required for Linear Solver= ',memuse2
!      Write(*,*) 'Total Gigabytes required for Newton Solver= ',
!     +               memuse1+memuse2
      ALLOCATE (A(Numnonzero), IA(Numcell+1))
      ALLOCATE (JA(Numnonzero))
      ALLOCATE (BB(Numcell), Hchange(Numcell))
      ALLOCATE (HCHOLD(Numcell),Wsave(Numcell))
      ALLOCATE (Dc(Numcell,6))
      A = 0.0D0
      Dc = 0.0D0
      BB = 0.0D0
      Hchange = 0.0D0
      Hchold = 1.0D0
      Wsave = 0.0D0 
      IA = 0
      JA = 0
!  
      ALLOCATE (Cvm1, Hvm1, Hvp1, Crm1, Hrm1, Hrp1, Ccm1)
      ALLOCATE (Hcm1, Hcp1, Ccc, Crr, Cvv, H)
      ALLOCATE (Hcoff, Rhss)
      ALLOCATE (W, Fhead, Fflux, Fheadsave, NJA)    
!
      W = 1.0D0
      Fhead = 0.0D0
      Fheadsave = 0.0D0
      Fflux = 0.0D0
! Order cells for jacobian and create CRS pointers
      CALL FILLINDEX(jj)
      NJA = IA(Numactive+1) - 1
      IF ( Linmeth.EQ.1 ) THEN
        CALL GMRES7AR(IN)
      ELSEIF ( Linmeth==2 ) THEN
        CALL XMD7AR(IN)
 !     ELSEIF ( Linmeth==3 ) THEN
!        CALL SAMG7AR(IN)
      END IF
!
!
C
! 
!
!--INITIALIZE SOLVER ARRAYS
      Cvm1 = 0.0D0
      Hvm1 = 0.0D0
      Crm1 = 0.0D0
      Hrm1 = 0.0D0
      Hrp1 = 0.0D0
      Ccm1 = 0.0D0
      Hcm1 = 0.0D0
      Hcp1 = 0.0D0
      Ccc = 0.0D0
      Crr = 0.0D0
      Cvv = 0.0D0
      H = 0.0D0
      Hcoff = 0.0D0
      Rhss = 0.0D0
!-------SAVE POINTERS FOR GRID AND RETURN
      CALL SGWF2NWT1PSV(Igrid)
!
      END SUBROUTINE GWF2NWT1AR
!     -----------------------------------------------------------------
!
      SUBROUTINE ORDERCELL()
! Order system for MODFLOW storage scheme by relating Row, Col, and Lay to
! Jacobian order
      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Iout, HNEW
      USE GWFBASMODULE, ONLY: HDRY
      USE GWFNWTMODULE
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER ic, ir, il, nrncnl, ij, jj, isum
!     -----------------------------------------------------------------
!
!ij points to a new row
!jj points to col. in sol. matrix
      ic = 1
      ir = 1
      il = 1
      ij = 1
      nrncnl = Nrow*Ncol*Nlay
      DO jj = 1, nrncnl
        IF ( ic.EQ.Ncol+1 ) THEN
          ic = 1
          ir = ir + 1
        ENDIF
        IF ( ir.EQ.Nrow+1 ) THEN
          ic = 1
          ir = 1
          il = il + 1
        ENDIF
        isum = 0
        IF ( Ibound(ic, ir, il).NE.0 ) THEN
          IF ( NCOL+NROW.LT.7 ) THEN
            IF ( IL.GT.1 ) isum = isum + abs(Ibound(ic, ir, il-1))
            IF ( IL.LT.NLAY ) isum = isum + abs(Ibound(ic, ir, il+1))
          END IF
          IF ( IR.GT.1 ) isum = isum + abs(Ibound(ic, ir-1, il))
          IF ( IC.GT.1 ) isum = isum + abs(Ibound(ic-1, ir, il))
          IF ( IR.LT.NROW ) isum = isum + abs(Ibound(ic, ir+1, il))
          IF ( IC.LT.NCOL ) isum = isum + abs(Ibound(ic+1, ir, il))
          IF ( isum.GT.0 ) THEN
            Diag(ij, 1) = il
            Diag(ij, 2) = ir
            Diag(ij, 3) = ic
            Icell(ic, ir, il) = ij
            ij = ij + 1
          ELSE
            WRITE(IOUT,*)
            WRITE(IOUT,*)'**Active cell surrounded by inactive cells**'
            WRITE(IOUT,*)'**Resetting cell to inactive**'
            WRITE(IOUT,*)'ROW=',ir,'COL=',ic,'LAY=',il
            WRITE(IOUT,*)
            Ibound(ic, ir, il) = 0
            HNEW(ic,ir,il) = HDRY
          END IF
        ENDIF
        ic = ic + 1
      ENDDO
      Numactive = ij - 1
      RETURN
      END SUBROUTINE ORDERCELL
!
!     -----------------------------------------------------------------
      SUBROUTINE FILLINDEX(jj)
! Fill CRS pointers for MODFLOW storage scheme by relating Row, Col, and Lay to
! Jacobian order
      USE GWFNWTMODULE
      USE GLOBAL, ONLY: IOUT,NCOL,NROW,NLAY,IBOUND
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER ic, ir, il, nrncnl, ij, jj
      INTEGER ILM, ILP, IRM, IRP, ICM, ICP, itemp
!     -----------------------------------------------------------------
!
!ij is the number of active cells (row in sol. vector)
!jj is the number of non-zero elements in the Jacobian
!IA() is the pointer to a new row in Jacobian (CRS)
!JA() points to the vector of unknowns for each entry in Jacobian (CRS)
      IA = 0
      JA = 0
      jj = 1
      DO ij = 1, Numactive
        il = Diag(ij, 1)
        ir = Diag(ij, 2)
        ic = Diag(ij, 3)
        IA(ij) = jj
! DIAGONAL FIRST
        JA(jj) = Icell(ic, ir, il)
        jj = jj + 1
        IF ( il.GT.1 ) THEN
          IF ( IBOUND(IC,IR,IL-1).NE.0 ) THEN
            JA(jj) = Icell(ic, ir, il-1)
            jj = jj + 1
          END IF
        ENDIF
        IF ( ir.GT.1 ) THEN
          IF ( IBOUND(IC,IR-1,IL).NE.0 ) THEN
            JA(jj) = Icell(ic, ir-1, il)
            jj = jj + 1
          END IF
        ENDIF
        IF ( ic.GT.1 ) THEN
          IF ( IBOUND(IC-1,IR,IL).NE.0 ) THEN
            JA(jj) = Icell(ic-1, ir, il)
            jj = jj + 1
          END IF
        ENDIF
        IF ( ic.LT.NCOL ) THEN
          IF ( IBOUND(IC+1,IR,IL).NE.0 ) THEN
            JA(jj) = Icell(ic+1, ir, il)
            jj = jj + 1
          END IF
        ENDIF
        IF ( ir.LT.NROW ) THEN
          IF ( IBOUND(IC,IR+1,IL).NE.0 ) THEN
            JA(jj) = Icell(ic, ir+1, il)
            jj = jj + 1
          END IF
        ENDIF
        IF ( il.LT.NLAY ) THEN
          IF ( IBOUND(IC,IR,IL+1).NE.0 ) THEN
            JA(jj) = Icell(ic, ir, il+1)
            jj = jj + 1
          END IF
        END IF
      ENDDO
      IA(Numactive+1) = jj
!      DO ij = 1, Numactive+1
!      write(iout,*)ia(ij)
!      end do
!      DO ij = 1, Numactive
!      DO jj = IA(ij),IA(ij+1)-1
!      write(iout,*)ja(jj)
!      end do 
!      end do
      RETURN     
      END SUBROUTINE FILLINDEX
!
!     -----------------------------------------------------------------
      SUBROUTINE COUNTACTIVE(jj)
      USE GWFNWTMODULE
      USE GLOBAL, ONLY: IOUT,NCOL,NROW,NLAY,IBOUND
      IMPLICIT NONE
!     ------------------------------------------------------------------
!     SPECIFICATIONS:
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     LOCAL VARIABLES
!     -----------------------------------------------------------------
      INTEGER ic, ir, il, nrncnl, ij, jj
      INTEGER ILM, ILP, IRM, IRP, ICM, ICP, itemp
!     -----------------------------------------------------------------
!
      jj = 0
      DO ij = 1, Numactive
        il = Diag(ij, 1)
        ir = Diag(ij, 2)
        ic = Diag(ij, 3)
! DIAGONAL FIRST
        jj = jj + 1
        IF ( il.GT.1 ) THEN
          IF ( IBOUND(IC,IR,IL-1).NE.0 ) THEN
            jj = jj + 1
          END IF
        ENDIF
        IF ( ir.GT.1 ) THEN
          IF ( IBOUND(IC,IR-1,IL).NE.0 ) THEN
            jj = jj + 1
          END IF
        ENDIF
        IF ( ic.GT.1 ) THEN
          IF ( IBOUND(IC-1,IR,IL).NE.0 ) THEN
            jj = jj + 1
          END IF
        ENDIF
        IF ( ic.LT.NCOL ) THEN
          IF ( IBOUND(IC+1,IR,IL).NE.0 ) THEN
            jj = jj + 1
          END IF
        ENDIF
        IF ( ir.LT.NROW ) THEN
          IF ( IBOUND(IC,IR+1,IL).NE.0 ) THEN
            jj = jj + 1
          END IF
        ENDIF
        IF ( il.LT.NLAY ) THEN
          IF ( IBOUND(IC,IR,IL+1).NE.0 ) THEN
            jj = jj + 1
          END IF
        END IF
      ENDDO
      RETURN     
      END SUBROUTINE COUNTACTIVE
!
      module unusedcode
!1     SUBROUTINE TEMPFILLUN. SET SCALERS FOR UNCONFINED FLOW
*      SUBROUTINE TEMPFILLUN(Ic, Ir, Il)
*      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Cv, Hnew, Cc, Cr, Ibound, Hcof,
*     +    Rhs, Botm, Lbotm,Iout
*      USE GWFNWTMODULE, ONLY:Cvm1,Hvp1,Hvm1,Crm1,Hrm1,Hrp1,Ccm1,Hcm1,
*     +                       Hcp1,Ccc,Crr,Cvv,H,Closezero,Icell,Hcoff,
*     +                       Rhss
*      USE GWFUPWMODULE, ONLY:Sn
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     ARGUMENTS
*!     ------------------------------------------------------------------
*      INTEGER Ic, Ir, Il, ij
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      INTEGER itemp
*      DOUBLE PRECISION THICK
*!     -----------------------------------------------------------------
*      Cvm1 = 0.0D0
*      Hvp1 = 0.0D0
*      Hvm1 = 0.0D0
*      Crm1 = 0.0D0
*      Hrm1 = 0.0D0
*      Hrp1 = 0.0D0
*      Ccm1 = 0.0D0
*      Hcm1 = 0.0D0
*      Hcp1 = 0.0D0
*      Ccc = 0.0D0
*      Crr = 0.0D0
*      Cvv = 0.0D0
*      H = Hnew(Ic, Ir, Il)
*      IF ( Ir.LT.Nrow ) THEN
*        IF ( IBOUND(IC,IR+1,IL).NE.0 ) THEN
*          Hrp1 = Hnew(Ic, Ir+1, Il)
*          IF ( Hrp1-H.GT.CLOSEZERO)THEN
*            THICK = dble(BOTM(IC,IR+1,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC,IR+1,LBOTM(IL)))
*            ij  = Icell(IC,IR+1,IL)
*            Ccc = Cc(Ic, Ir, Il)*THICK*Sn(ij)
*          ELSE
*            THICK = dble(BOTM(IC,IR,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC,IR,LBOTM(IL)))
*            ij  = Icell(IC,IR,IL)
*            Ccc = Cc(Ic, Ir, Il)*THICK*Sn(ij)
*          END IF
*        END IF
*      ENDIF
*      IF ( Ic.LT.Ncol ) THEN
*        IF ( IBOUND(IC+1,IR,IL).NE.0 ) THEN
*          Hcp1 = Hnew(Ic+1, Ir, Il)
*          IF ( Hcp1-H.GT.CLOSEZERO )THEN
*            THICK = dble(BOTM(IC+1,IR,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC+1,IR,LBOTM(IL)))
*            ij = Icell(IC+1,IR,IL)
*            Crr = Cr(Ic, Ir, Il)*THICK*Sn(ij)
*          ELSE 
*            THICK = dble(BOTM(IC,IR,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC,IR,LBOTM(IL)))
*            ij = Icell(IC,IR,IL)
*            Crr = Cr(Ic, Ir, Il)*THICK*Sn(ij)
*          END IF
*        END IF
*      ENDIF
* ! Need to correct the following calculations for case when CV is head dependent.
*      IF ( Il.LT.Nlay ) THEN
*        IF ( IBOUND(IC,IR,IL+1).NE.0 ) THEN
*          Hvp1 = Hnew(Ic, Ir, Il+1)
* ! Set vertical correction for perched conditions
*!          IF ( Hvp1.LT.BOTM(IC,IR-1,LBOTM(IL) ) 
*!    +         Hvp1 = BOTM(IC,IR-1,LBOTM(IL)
*          Cvv = Cv(Ic, Ir, Il)
*        END IF
*      ENDIF
*      IF ( Il.GT.1 ) THEN
*        IF ( IBOUND(IC,IR,IL-1).NE.0 ) THEN
*          Hvm1 = Hnew(Ic, Ir, Il-1)
*          Cvm1 = Cv(Ic, Ir, Il-1)
*        END IF
*      ENDIF
*      IF ( Ir.GT.1 ) THEN
*        IF ( IBOUND(IC,IR-1,IL).NE.0 ) THEN
*          Hrm1 = Hnew(Ic, Ir-1, Il)
*          IF ( Hrm1-H.GT.CLOSEZERO)THEN
*            THICK = dble(BOTM(IC,IR-1,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC,IR-1,LBOTM(IL)))
*            ij = Icell(IC,IR-1,IL)
*            Ccm1 = Cc(Ic, Ir-1, Il)*THICK*Sn(ij)
*          ELSE
*            THICK = dble(BOTM(IC,IR,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC,IR,LBOTM(IL)))
*            ij = Icell(IC,IR,IL)
*            Ccm1 = Cc(Ic, Ir-1, Il)*THICK*Sn(ij)
*          END IF
*        END IF
*      ENDIF
*      IF ( Ic.GT.1 ) THEN
*        IF ( IBOUND(IC-1,IR,IL).NE.0 ) THEN
*          Hcm1 = Hnew(Ic-1, Ir, Il)
*          IF ( Hcm1-H.GT.CLOSEZERO )THEN
*            THICK = dble(BOTM(IC-1,IR,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC-1,IR,LBOTM(IL)))
*            ij = Icell(IC-1,IR,IL)
*            Crm1 = Cr(Ic-1, Ir, Il)*THICK*Sn(ij)
*          ELSE 
*            THICK = dble(BOTM(IC,IR,LBOTM(IL)-1)) - 
*     +              dble(BOTM(IC,IR,LBOTM(IL)))
*            ij = Icell(IC,IR,IL)
*            Crm1 = Cr(Ic-1, Ir, Il)*THICK*Sn(ij)
*          END IF
*        END IF
*      ENDIF
*      Hcoff = Hcof(Ic, Ir, Il)
*      Rhss = Rhs(Ic, Ir, Il)
*      END SUBROUTINE TEMPFILLUN
!
!     -----------------------------------------------------------------
!
!     SUBROUTINE TEMPFILLCON. SET SCALERS FOR CONFINED FLOW
*      SUBROUTINE TEMPFILLCON(Ic, Ir, Il)
*      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Cv, Hnew, Cc, Cr, Ibound, Hcof,
*     +    Rhs, Botm, Lbotm,Iout
*      USE GWFNWTMODULE, ONLY:Cvm1,Hvp1,Hvm1,Crm1,Hrm1,Hrp1,Ccm1,Hcm1,
*     +                       Hcp1,Ccc,Crr,Cvv,H,Closezero,Icell,Hcoff,
*     +                       Rhss
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     ARGUMENTS
*!     ------------------------------------------------------------------
*      INTEGER Ic, Ir, Il, ij
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      INTEGER itemp
*      DOUBLE PRECISION THICK
*!     -----------------------------------------------------------------
*      Cvm1 = 0.0D0
*      Hvp1 = 0.0D0
*      Hvm1 = 0.0D0
*      Crm1 = 0.0D0
*      Hrm1 = 0.0D0
*      Hrp1 = 0.0D0
*      Ccm1 = 0.0D0
*      Hcm1 = 0.0D0
*      Hcp1 = 0.0D0
*      Ccc = 0.0D0
*      Crr = 0.0D0
*      Cvv = 0.0D0
*      H = Hnew(Ic, Ir, Il)
*      IF ( Ir.LT.Nrow ) THEN
*        IF ( IBOUND(IC,IR+1,IL).NE.0 ) THEN
*           Hrp1 = Hnew(Ic, Ir+1, Il)
*           Ccc = Cc(Ic, Ir, Il)
*        END IF
*      ENDIF
*      IF ( Ic.LT.Ncol ) THEN
*        IF ( IBOUND(IC+1,IR,IL).NE.0 ) THEN
*          Hcp1 = Hnew(Ic+1, Ir, Il)
*          Crr = Cr(Ic, Ir, Il)
*        END IF
*      ENDIF
* ! Need to correct the following calculations for case when CV is head dependent.
*      IF ( Il.LT.Nlay ) THEN
*        IF ( IBOUND(IC,IR,IL+1).NE.0 ) THEN
*          Hvp1 = Hnew(Ic, Ir, Il+1)
*          Cvv = Cv(Ic, Ir, Il)
*        END IF
*      ENDIF
*      IF ( Il.GT.1 ) THEN
*        IF ( IBOUND(IC,IR,IL-1).NE.0 ) THEN
*          Hvm1 = Hnew(Ic, Ir, Il-1)
*          Cvm1 = Cv(Ic, Ir, Il-1)
*        END IF
*      ENDIF
*      IF ( Ir.GT.1 ) THEN
*        IF ( IBOUND(IC,IR-1,IL).NE.0 ) THEN
*          Hrm1 = Hnew(Ic, Ir-1, Il)
*          Ccm1 = Cc(Ic, Ir-1, Il)
*        END IF
*      ENDIF
*      IF ( Ic.GT.1 ) THEN
*        IF ( IBOUND(IC-1,IR,IL).NE.0 ) THEN
*          Hcm1 = Hnew(Ic-1, Ir, Il)
*          Crm1 = Cr(Ic-1, Ir, Il)
*        END IF
*      ENDIF
*      Hcoff = Hcof(Ic, Ir, Il)
*      Rhss = Rhs(Ic, Ir, Il)
*      END SUBROUTINE TEMPFILLCON
!
!     -----------------------------------------------------------------
!
*      FUNCTION DHORIZ(Hup, Ttop, Bbot, il)
*! RETURNS DERIVATIVE OF HORIZONTAL CONDUCTANCE BASED ON SMOOTH FUNCTION
*! FUNCTION IS CALCULATED IN UPW PACKAGE IN SUBROUTINE SAT_THICK
*      USE GWFNWTMODULE
*      USE GWFUPWMODULE, ONLY: LAYTYPUPW
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     ------------------------------------------------------------------
*!     ------------------------------------------------------------------
*!     ARGUMENTS
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION Hup, Ttop, Bbot
*      DOUBLE PRECISION DHORIZ
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION factor, x, s, v, cof1, cof2, EPS, ACOF, Y
*      DOUBLE PRECISION EPSQD, z    
*      INTEGER il
*!     -----------------------------------------------------------------
*      DHORIZ = 0.0D0
*      IF ( LAYTYPUPW(il).LE.0 ) RETURN
*C-------STRAIGHT LINE WITH PARABOLIC SMOOTHING
*      EPS = Thickfact
*      ACOF = 1.0 / (1.0 - EPS)
*      x = (Hup-bbot)/(TTOP-BBOT)
*      IF ( x.LT.1.0d-9 )  x = 1.0d-9
*      IF(X.LT.EPS)THEN
*        Y = ACOF * X / (EPS*(Ttop-Bbot))
*      ELSEIF(X.LT.1.0D0-EPS)THEN
*        Y = ACOF /(Ttop-Bbot)
*      ELSEIF(X.LT.1.0D0)THEN
*        X = 1.0 - X
*        Y = - ACOF * x / (EPS * (Ttop - Bbot))
*        Y = 1.0-Y
*      ELSE
*        Y = 0.0
*      ENDIF
*      factor = Y
*      DHORIZ = factor     
*      END FUNCTION DHORIZ
!
!     -----------------------------------------------------------------
*      DOUBLE PRECISION FUNCTION DVERT(Hh,Ttop,Bbot)
*      USE GWFNWTMODULE
*      IMPLICIT NONE
*!     ******************************************************************
*!     COMPUTE THE DERIVATIVE OF THE VERTICAL BRANCH CONDUCTANCE BETWEEN
*!     A LAYER AND THE NEXT LOWER LAYER FROM VERTICAL HYDRAULIC 
*!     CONDUCTIVITY.
*!     ******************************************************************
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     ------------------------------------------------------------------
*!     ARGUMENTS
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION Hh,Ttop,Bbot
*      INTEGER J, I, K, iltyp
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION factor, thick, zero, x
*!     ------------------------------------------------------------------
*!
*      DVERT = 0.0D0
*      RETURN   ! CV is constant as of now.
*      zero = 0.0D0
*      factor = 1.0D0
*      thick = Thickfact*(ttop - bbot)
*      x = Hh-bbot  
*      factor = 1.0D0/thick
*      IF ( x.LE.0.0D0 ) factor = 0.0D0
*      IF ( x.GT.thick ) factor = 0.0D0
*      DVERT = factor
*      END FUNCTION DVERT
!
!     -----------------------------------------------------------------
!
!     -----------------------------------------------------------------
*      SUBROUTINE GWF2NWT1FM(Kkiter, ICNVG, KSTP, KPER, Maxiter, 
*     +                      Iunitchd, Igrid)
*! Builds and Solves Jacobian
*! Calls various unstructured linear solvers to solve Jacobian
*      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Hcof, Rhs, Iout,botm,
*     +                 LBOTM, HOLD, HNEW, DELR, DELC, ISSFLG
*      USE GWFBASMODULE, ONLY:TOTIM, HNOFLO
*      USE GWFNWTMODULE
*      USE XMDMODULE
*      USE GMRESMODULE
*!      USE SAMGMODULE
*      USE ilupc_mod
*      USE machine_constants
*      use SimPHMFModule, only: ustop
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     ------------------------------------------------------------------
*      INTEGER, EXTERNAL :: IACTIVE
*      INTEGER, EXTERNAL :: IGETUNIT
*      DOUBLE PRECISION, EXTERNAL :: GW_func, RMS_func
*      EXTERNAL SGWF2NWT1PNT, TEMPFILLUN, TEMPFILLCON
*!     ------------------------------------------------------------------
*!     ARGUMENTS
*!     ------------------------------------------------------------------
*      INTEGER Iss, Igrid, Kkiter, Icnvg, Maxiter, KSTP, KPER, Iunitchd
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION h2
*      DOUBLE PRECISION FHEADTEMP,R_norm_gmres
*      DOUBLE PRECISION r_norm, fheadsave2, STPNUM
*      INTEGER ic, ir, il, LICNVG,ITER, ippn
*      INTEGER ij, jj, ichld, irhld, ilhld, IUNSTR, itertot
*      INTEGER n_iter, z, n, icfld, irfld, ilfld, ITP, ITER1U
*!      CHARACTER ierr
*      SAVE FHEADTEMP,itertot
*      DOUBLE PRECISION TINY,SCALE
*
*      INTEGER matrix,ncyc_done,ierr_samg,control,iout_samg,i,j,new
*      DOUBLE PRECISION res_out
*!     -----------------------------------------------------------------
*!
*!
*!1------SET POINTERS FOR THE CURRENT GRID.
*      CALL SGWF2NWT1PNT(Igrid)
*      TINY=1.0D-30
*      fheadsave2 = 0.0D0
*      ic = 1
*      ir = 1
*      il = 1
*      ij = 0
*      ichld = 1
*      irhld = 1
*      ilhld = 1
*      icfld = 1
*      irfld = 1
*      ilfld = 1
*      Icnvg = 0
*      iierr = 0
*      ippn = 0
*      r_norm = 0.0D0
*      n_iter = 0
*      ISS=ISSFLG(KPER)
*      CALL HEAD_SAVE()
*! Save head as previous iteration (Hiter)
*!     SOLVE FOR GROUNDWATER HEAD IN 3D USING NEWTON/PICARD
*      IF ( Kkiter==1 ) THEN
*        itertot = 0
*        II = 0
*        RMS2 = 0.0
*        IF ( kkiter+kper.EQ.2) Fheadsave = 2.0*Tol
*        Fhead = 0.0D0
*!        IF ( Iunitchd.GT.0 ) THEN
*!          CALL ORDERCELL()
*!          CALL FILLINDEX(jj)
*!        END IF
*      END IF
*      IF ( II.EQ.0 ) RMSAVE = RMS1
*      RMS2 = RMS1
*      RMS1 = RMS_func(icfld,irfld,ilfld)
*      Icnvg = 0
*      IF ( RMS1.GT.FTOL .OR. ABS(Fheadsave).GT.Tol .OR. 
*     +                           kkiter.LT.2 ) THEN
*        Ibt = 1
*        IF ( BTRACK.EQ.0 .OR. II.GE.Numtrack ) Ibt = 0
*        IF ( RMS1.LT.Btol*rmsave .OR. Kkiter.EQ.1 ) Ibt = 0
*        IF ( II.GT.0 .AND. RMS1.GT.RMS2 ) Ibt = 0
*        IF ( Ibt.EQ.0 ) THEN
*          II = 0
*          jj = 1
*!          IF ( Nonmeth==1 ) THEN
*            CALL Jacobian(kkiter,kper,kstp)
*!          ELSE
*!            CALL Picard(ij,jj,kkiter)
*!          END IF
*! Calls for linear solver
*          n = Numactive
*          IF ( Linmeth.EQ.1 ) THEN
*! M_save_dir is size of Krylov space. Same as msdr:
*            SELECT CASE (ilu_method)
*            CASE (1)
*            CALL ilut (n,A,JA,IA,Lev_fill,Drop_tol,Alu,Jlu,Ju,Iierr)
*            CASE (2)
*            CALL iluk (n,A,JA,IA,lev_fill,Alu,Jlu,Ju,Iierr)
*            END SELECT
*            IF(Iierr /= 0) THEN
*            WRITE(Iout,*) 'Error in Preconditioning: ', Iierr
*            CALL USTOP('  ')
*            END IF
*            CALL gmres(n,Msdr,BB,Hchange,Stop_tol_gmres,Maxitr_gmres,A,
*     +             JA,IA,Alu,Jlu,Ju,iierr,n_iter,R_norm_gmres)
*            IF(Iierr == 1) THEN
*!            WRITE(Iout,*) 'Linear solver failed to converge: '
*!     +                  , Iierr,n_iter, R_norm_gmres
*!            CALL USTOP('  ')
*            ELSEIF(Iierr /= 0) THEN  
*              WRITE(Iout,*) 'Error in gmres: ', Iierr,n_iter, 
*     +                       R_norm_gmres
*              CALL USTOP('  ')
*            END IF
*          ELSEIF(LINMETH.EQ.2)THEN
*C
*C---------CALL XMD SOLVER     
*c
*            IF (IDROPTOL.EQ.0 .or.  kkiter.gt.1) THEN
*              call xmdnfctr(a, bb, ia, ja, nja, numactive, ierr)
*            ELSE
*              call xmdprecd(a, bb, epsrn, ia, ja, nja, numactive,
*     [                      level, ierr)
*
*cmi
*C
*            ENDIF                                                     
*c  -----------------
*c         solve matrix
*          iter = Mxiterxmd
*          call xmdsolv(a, bb, hchange, hclosexmd, rrctol, ia, ja, nja,
*     [                 numactive, north, iter, iacl, ierr)
*            n_iter = iter
*C     
*!          ELSE IF ( Linmeth.EQ.3 ) THEN !SAMG solver
*
*c         INPUT DATA
*c            n=number of cells
*c            nja=number of non-zero entries in the linear system
*c            ia,ja,a matrix in CRS format
*c            BB= right hand side of the linear system
*c            Hchange= first guess / solution of the linear system
*c            matrix=22 if not symmetric else matrix=12
*!             matrix=22 
*c            Stop_tol_samg = closure criterion L2-norm 
*c               if eps<0 :: absolute residual 
*c               else     :: relative reduction compared to first guess
*c            Maxitr_samg=maximal number of iterations
*c            control=1 no setup reuse
*c            control=2 automatic setup reuse
*!             control=1
*c            iout  Controls print output related to SAMG�s solution phase.
*c                <0 No printout, except for warnings and errors24.
*c                =0 Minimal output on results and timings.
*c                >0 Additional print output specified by the individual digits:
*c                   iout[1:1]
*c                      1 Table of input data and work statistics.
*c                      2 In addition: history of cycling process.
*c                      3 Extended history: including all levels, full smoothing steps.
*c                      4 Extended history: including all levels, partial smoothing steps.
*c                    iout[2:2]
*c                      0 No action.
*c                      1 Display most relevant hidden parameters.
*c                      2 Display all hidden parameters.
*c                      3 Display all hidden parameters in a single list.
*!              iout_samg=2
*
*c
*c         OUTPUT DATA
*c            res_out=output residual L2-Norm (Norm can be changed)
*c            ncyc_done=iterations performed
*c            ierr=error flag         
*             
*
*
*c !Diagonals have to be positive!
*!             do i=1,n
*!                if (a(ia(i)).lt.0) then
*!                   do j=ia(i),ia(i+1)-1
*!                      a(j)=-a(j)
*!                   enddo
*!                   bb(i)=-bb(i)
*!                endif
*!             enddo
*!
*!100          call samg_interface(n,nja,ia,ja,a,BB,Hchange, 
*!     +              matrix,res_out,ncyc_done,ierr_samg,          
*!     +              samg_approach,Stop_tol_samg,
*!     +              Maxitr_samg,control,iout_samg)
*!
*!!    DEBUG OPTION
*!             write(6,*) res_out,ierr_samg,ncyc_done
*!!             if (ierr_samg.eq.0) goto 100
*!
*!             if (ierr_samg.gt.0) then
*!              WRITE(Iout,*) 'Error in SAMG - err,cyc,res:', ierr_samg,
*!     +                        ncyc_done,res_out
*!              CALL USTOP('  ')
*!             endif
*                
*          END IF
*C
*C--Update heads.
*          CALL GWF2NWT1UPH2(ichld,irhld,ilhld,ISS,Kkiter)
*          Fheadsave = Fhead
*          ippn = 0
*          IF ( RMS1.LT.FTOL .AND. ABS(Fheadsave).LT.Tol ) THEN
*            Icnvg = 1
*            ippn = 1
*          END IF
*          Itreal = Itreal + 1
*        ELSE
*          n_iter = 0
*          CALL Back_track(ichld,irhld,ilhld,ISS)
*          II = II + 1
*        END IF
* !       Itreal = Itreal + 1
*      ELSE
*        Icnvg = 1
* !     open(999,file='Init_head.dat')
* !     if(kper==1)then
* !     Do il = 1, Nlay
* !       Do ir = 1, Nrow
* !         write(999,888)(HNEW(ic,ir,il),ic = 1, Ncol)
* !       end do
* !     end do
* !     end if
*      END IF
* ! 888 format(256E20.10)
*!
*!  Calculate maximum head change and residuals
*!  Write head and flux residuals
*!  Write iteration header
*      IF ( IPRNWT.GT.0 ) THEN
*        IF ( MOD(Kkiter,10).eq.1 .AND. Icnvg.EQ.0 ) THEN
*          IF ( Btrack.GT.0 ) THEN
*            WRITE(IOUT,111)
*          ELSE
*            WRITE(IOUT,112)
*          END IF
*        END IF
*      END IF
*  111 FORMAT (1X,'                                  ',
*     +           '                Max.-Head-Change',
*     +           '                        ',
*     +           'Max.-Flux-Residual',/
*     +          '  Residual-Control   Outer-Iter.   ',
*     +          'Inner-Iter.    ',
*     +          ' Column Row Layer ',   
*     +          '   Max.-Head-Change   ', 
*     +          ' Column Row Layer ',
*     +          '   Max.-Flux-Residual',
*     +          '            L2-New               L2-Old           ',
*     +          'Solver-Max-Delh' )
*  112 FORMAT (1X,'                                ',
*     +           '                  Max.-Head-Change',
*     +           '                         ',
*     +           'Max.-Flux-Residual',/
*     +          '  Residual-Control   Outer-Iter.   ',
*     +          'Inner-Iter.    ',
*     +          ' Column Row Layer ',   
*     +          '   Maximum-Head-Change ', 
*     +          ' Column Row Layer ',
*     +          '   Maximum-Flux-Residual          L2-NORM ')
*      itertot = itertot + n_iter
*      IF ( IPRNWT.GT.0 ) THEN
*        IF ( Icnvg.EQ.0 .OR. ippn.EQ.1) THEN
*          IF ( Btrack.GT.0 ) THEN 
*                      WRITE (Iout, 9001) II,itreal,n_iter,
*     +                 ichld,irhld,ilhld,fhead,icfld,irfld,ilfld,
*     +                 Fflux,RMS1,RMS2,FHEADSAVE
*          ELSE
*            WRITE (Iout, 9002) II,itreal,n_iter,ichld,irhld,ilhld,fhead,
*     +                   icfld,irfld,ilfld,fflux,RMS1
* !    +      HNEW(ichld,irhld,ilhld),BOTM(ichld,irhld,ilhld-1),
* !    +      BOTM(ichld,irhld,ilhld)
*          END IF
*        END IF
*      END IF
*      IF ( Icnvg.GT.0 .OR. itreal.GE.Maxiter) THEN
*        WRITE (Iout,9003) itreal, itertot
*      END IF
* 9001 FORMAT (5X,I6,12X,I6,6X,I6,8X,I6,1x,I4,3X,I3,3X,E20.10,
*     +        2x,I6,1x,I4,3X,I2,1X,4(2X,E20.10))
* 9002 FORMAT (5X,I6,12X,I6,6X,I6,8X,I6,1x,I4,3X,I3,3X,E20.10,
*     +        2x,I6,1x,I4,3X,I2,3X,2(2X,E20.10))
* 9003 FORMAT (/4x,'------------------------------------------------',/
*     +        7x,'NWT REQUIRED     ',i8,' OUTER ITERATIONS ',/
*     +        7x,'AND A TOTAL OF   ',i8,' INNER ITERATIONS.',/
*     +        4x,'------------------------------------------------')
*      END SUBROUTINE GWF2NWT1FM
!
!
!     -----------------------------------------------------------------
!     Set the head in dewatered cells to Hdry.
*      SUBROUTINE GWF2NWT1BD()
*      USE GLOBAL, ONLY:Hnew, Nrow, Ncol, Nlay, BOTM, LBOTM, IBOUND,
*     +                 LAYHDT, IOUT  
*      USE GWFBASMODULE,ONLY:HDRY
*      USE GWFUPWMODULE,ONLY:IPHDRY
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      INTEGER ic, ir, il
*!     -----------------------------------------------------------------
*!      DO ir = 1, Nrow
*!        write(iout,222)( Hnew(ic, ir, 1), ic = 1, Ncol)
*!      end do
*! 222  format(113e20.10)
*      CALL Head_save()   !From Scott B. 9/7/2013
*      DO il = 1, Nlay
*        IF ( LAYHDT(il).GT.0 ) THEN
*          DO ir = 1, Nrow
*            DO ic = 1, Ncol
*              IF ( IBOUND(ic,ir,il).GT.0 .AND. IPHDRY.GT.0 ) THEN
*                IF ( Hnew(ic, ir, il)-dble(BOTM(ic,ir,LBOTM(il)))
*     +                                                  .LT.2.0e-3 )
*     +               Hnew(ic, ir, il) = dble(Hdry)
*              END IF
*            ENDDO
*          ENDDO
*        END IF
*      ENDDO
*!      do ir=1,nrow
*!      write(iout,101)(Hnew(ic, ir, 1),ic=1,ncol)
*!      end do
*!  101 format(80E20.10) 
*      END SUBROUTINE GWF2NWT1BD
!
!
!     -----------------------------------------------------------------
!
*      SUBROUTINE Back_track(ichld,irhld,ilhld,ISS)
*      USE GLOBAL, ONLY:Ibound, Hnew, Lbotm, Botm, Iout, NLAY
*      USE GWFUPWMODULE, ONLY:LAYTYPUPW
*      USE GWFNWTMODULE
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*      DOUBLE PRECISION, EXTERNAL :: Sum_sat
*!     -----------------------------------------------------------------
*!
*!     -----------------------------------------------------------------
*!     ARGUMENTS
*!     -----------------------------------------------------------------
*      INTEGER ichld,irhld,ilhld,ISS
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      INTEGER ic, ir, il, jj, i, ibotlay
*      DOUBLE PRECISION hsave, sum
*!     -----------------------------------------------------------------
*!  APPLY BACKTRACKING
*!  
*      Fhead = 0.0D0    
*      DO jj = 1, Numactive
*        il = Diag(jj, 1)
*        ir = Diag(jj, 2)
*        ic = Diag(jj, 3)
*        Hiter(ic, ir, il) = Hiter(ic, ir, il) - Hchange(jj)
*        Hchange(jj) = Breduc*Hchange(jj)
*        HNEW(ic, ir, il) = HITER(ic, ir, il) + Hchange(jj)
*        IF ( IBOTAV.GT.0 .AND. LAYTYPUPW(il).GT.0 ) THEN
*          ibotlay = il
*          DO i = il + 1, NLAY - 1
*            IF ( IBOUND(ic,ir,i).GT.0 ) ibotlay = ibotlay + 1
*          END DO
*          IF ( il.EQ.NLAY ) THEN
*            IF (Hnew(ic, ir, il).LT.dble(Botm(ic, ir, Lbotm(il))) ) THEN
*              IF ( Hiter(ic, ir, il).LT.dble(BOTM(ic, ir, Lbotm(il))) )
*     +          Hiter(ic, ir, il) = BOTM(ic, ir, Lbotm(il)) + 
*     +                              1.0e-6
*              sum = 0.0D0
*              sum = Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat
*     +              (Sum_sat(sum,ic,ir,il),ic-1,ir,il),ic+1,ir,il),
*     +               ic,ir-1,il),ic,ir+1,il),ic,ir,il-1),ic,ir,il+1)
*              IF ( sum .LT.1.0e-7 ) THEN
*                hsave = Hnew(ic, ir, il)
*                Hnew(ic, ir, il) = (Hiter(ic, ir, il)+
*     +                              dble(BOTM(ic, ir, Lbotm(il))))/2.0d0
*                Hchange(jj) = Hnew(ic, ir, il) - hsave
*              END IF
*            END IF
*          ELSEIF ( IBOUND(ic,ir,ibotlay+1).EQ.0 ) THEN
*            IF (Hnew(ic, ir, il).LT.dble(Botm(ic, ir, Lbotm(ibotlay))) )
*     +                                                        THEN
*              IF ( Hiter(ic, ir, il).LT.
*     +                              dble(BOTM(ic, ir, Lbotm(ibotlay))) )
*     +          Hiter(ic, ir, il) = dble(BOTM(ic, ir, Lbotm(ibotlay))) +
*     +                              1.0e-6
*              sum = 0.0D0
*              sum = Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat
*     +              (Sum_sat(sum,ic,ir,il),ic-1,ir,il),ic+1,ir,il),
*     +               ic,ir-1,il),ic,ir+1,il),ic,ir,il-1),ic,ir,il+1)
*              IF ( sum .LT.1.0e-7 ) THEN
*                hsave = Hnew(ic, ir, il)
*                Hnew(ic, ir, il) = (Hiter(ic, ir, il)+
*     +                     dble(BOTM(ic, ir, Lbotm(ibotlay))))/2.0d0
*                Hchange(jj) = Hnew(ic, ir, il) - hsave
*              END IF
*            END IF
*          ENDIF
*        END IF
*        IF ( ABS(Hchange(jj)).GT.ABS(fhead) ) THEN
*          Fhead = Hchange(jj)
*          ichld = ic
*          irhld = ir
*          ilhld = il
*        ENDIF
*      END DO
*      END SUBROUTINE Back_track
!
!
!     -----------------------------------------------------------------
*      SUBROUTINE GWF2NWT1UPH2(ichld,irhld,ilhld,ISS,Kkiter)
*!  Update heads and apply relaxation (delta-bar-delta)
*      USE GLOBAL, ONLY:Ibound, Hnew, Lbotm, Botm, Iout, NCOL,NROW,NLAY
*      USE GWFUPWMODULE, ONLY:LAYTYPUPW
*      USE GWFNWTMODULE
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     -----------------------------------------------------------------
*      INTRINSIC ABS
*      DOUBLE PRECISION, EXTERNAL :: GW_func, Sum_sat
*!     -----------------------------------------------------------------
*!     ARGUMENTS
*!     -----------------------------------------------------------------
*      INTEGER Kkiter,Icnvg,ichld,irhld,ilhld,ISS
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION s, wstar, ww, change, hsave, sum
*      INTEGER ic, ir, il, jj, ibotlay, i
*!     -----------------------------------------------------------------
*!
*      Fhead = 0.0D0
*      DO jj = 1, Numactive
*        il = Diag(jj, 1)
*        ir = Diag(jj, 2)
*        ic = Diag(jj, 3)
*        Hchange(jj) = Hchange(jj) - Hnew(ic,ir,il)
*        IF ( kkiter.EQ.1 )THEN
*          Wsave(jj) = 1.0D0
*          Hchold(jj) = Hchange(jj)
*        END IF            
*        ww = Wsave(jj)
*        IF ( Hchold(jj)*Hchange(jj).LT.0.0D0 ) THEN
*          ww = Theta*Wsave(jj)
*        ELSE
*          ww = Wsave(jj) + akappa
*        END IF
*        IF ( ww.GT.1.0d0 ) ww = 1.0d0
*        Hchold(jj) = (1-gamma) * Hchange(jj) + gamma * Hchold(jj)
*        Wsave(jj) = ww
*        Hchange(jj) = Hchange(jj) * ww + amomentum * Hchold(jj)
*        Hnew(ic, ir, il) = Hiter(ic, ir, il) + Hchange(jj)
*        IF ( IBOTAV.GT.0 .AND. LAYTYPUPW(il).GT.0 ) THEN
*          ibotlay = il
*          DO i = il + 1, NLAY - 1
*            IF ( IBOUND(ic,ir,i).GT.0 ) ibotlay = ibotlay + 1
*          END DO
*          IF ( il.EQ.NLAY ) THEN
*            IF (Hnew(ic, ir, il).LT.dble(Botm(ic, ir, Lbotm(il))) ) THEN
*              IF ( Hiter(ic, ir, il).LT.dble(BOTM(ic, ir, Lbotm(il))) )
*     +          Hiter(ic, ir, il) = BOTM(ic, ir, Lbotm(il)) + 
*     +                              1.0e-6
*              sum = 0.0D0
*              sum = Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat
*     +              (Sum_sat(sum,ic,ir,il),ic-1,ir,il),ic+1,ir,il),
*     +               ic,ir-1,il),ic,ir+1,il),ic,ir,il-1),ic,ir,il+1)
*              IF ( sum .LT.1.0e-7 ) THEN
*                hsave = Hnew(ic, ir, il)
*                Hnew(ic, ir, il) = (Hiter(ic, ir, il)+
*     +                              dble(BOTM(ic, ir, Lbotm(il))))/2.0d0
*                Hchange(jj) = Hnew(ic, ir, il) - hsave
*              END IF
*            END IF
*          ELSEIF ( IBOUND(ic,ir,ibotlay+1).EQ.0 ) THEN
*            IF (Hnew(ic, ir, il).LT.dble(Botm(ic, ir, Lbotm(ibotlay))) )
*     +                                                              THEN
*              IF ( Hiter(ic, ir, il).LT.
*     +                              dble(BOTM(ic, ir, Lbotm(ibotlay))) )
*     +          Hiter(ic, ir, il) = dble(BOTM(ic, ir, Lbotm(ibotlay))) +
*     +                              1.0e-6
*              sum = 0.0D0
*              sum = Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat(Sum_sat
*     +              (Sum_sat(sum,ic,ir,il),ic-1,ir,il),ic+1,ir,il),
*     +               ic,ir-1,il),ic,ir+1,il),ic,ir,il-1),ic,ir,il+1)
*              IF ( sum .LT.1.0e-7 ) THEN
*                hsave = Hnew(ic, ir, il)
*                Hnew(ic, ir, il) = (Hiter(ic, ir, il)+
*     +                    dble(BOTM(ic, ir, Lbotm(ibotlay))))/2.0d0
*                Hchange(jj) = Hnew(ic, ir, il) - hsave
*              END IF
*            END IF
*          ENDIF
*        END IF
*        IF ( ABS(Hchange(jj)).GT.ABS(fhead) ) THEN
*          fhead = Hchange(jj)
*          ichld = ic
*          irhld = ir
*          ilhld = il
*        ENDIF
*      ENDDO
*      END SUBROUTINE GWF2NWT1UPH2
!
!
!     -----------------------------------------------------------------
!     Save previous iteration heads.
*      SUBROUTINE Head_save()
*      USE GLOBAL, ONLY:Hnew, Nrow, Ncol, Nlay, Iout, BOTM
*      USE GWFBASMODULE,ONLY:HDRY
*      USE GWFNWTMODULE
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      INTEGER ic, ir, il
*!     -----------------------------------------------------------------
*      DO il = 1, Nlay
*        DO ir = 1, Nrow
*          DO ic = 1, Ncol
*            Hiter(ic, ir, il) = Hnew(ic, ir, il)
*          ENDDO
*        ENDDO
*      ENDDO
*      END SUBROUTINE HEAD_SAVE
!
!     -----------------------------------------------------------------
!     Return value of groundwater flow equation
*      DOUBLE PRECISION FUNCTION GW_func(Ic, Ir, Il)
*      USE GWFNWTMODULE
*      USE GLOBAL,      ONLY:iout
*      USE GWFBASMODULE, ONLY:HNOFLO
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     -----------------------------------------------------------------
*!     ARGUMENTS
*!     -----------------------------------------------------------------
*      INTEGER Ic, Ir, Il, Iunitlak
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION term1, term2, term3
*!     -----------------------------------------------------------------   
*      GW_func = 0.0D0
*      IF ( H==HNOFLO ) RETURN    
*      term1 = Cvm1*Hvm1 + Ccm1*Hrm1 + Crm1*Hcm1
*      term2 = (-Cvm1-Ccm1-Crm1-Crr-Ccc-Cvv+Hcoff)*H
*      term3 = Crr*Hcp1 + Ccc*Hrp1 + Cvv*Hvp1 - Rhss
*      GW_func = term1 + term2 + term3
*!      if(ic==166.and.ir==159.and.il==1)then
*!      write(iout,222)ic,ir,il,cvm1*(Hvm1-h),ccm1*(hrm1-h),crm1*(hcm1-h),
*!     +cvv*(hvp1-h),ccc*(hrp1-h),crr*(hcp1-h),rhss,gw_func
*!      end if
*!  222 format(3i5,8e20.10)
*      !if(ibound(ic,ir,il)==-1)then
*      !if(ibound(ic,ir-1,il)==1)sum=sum+ccm1*(hrm1-h)
*      !if(ibound(ic-1,ir,il)==1)sum=sum+crm1*(hcm1-h)
*      !if(ibound(ic+1,ir,il)==1)sum=sum+crr*(hcp1-h)
*      !if(ibound(ic,ir+1,il)==1)sum=sum+ccc*(hrp1-h)
*      !write(iout,*)'sum=',ic,ir,il,sum
*      !end if
*      END FUNCTION GW_func
!
!
!
!     -----------------------------------------------------------------
!     Return value of L2-Norm of GW equation, max flux
!     and head residuals
*      DOUBLE PRECISION FUNCTION RMS_func(icfld,irfld,ilfld)
*      USE GWFNWTMODULE, ONLY: Fflux,Numactive,Diag
*      USE GWFUPWMODULE, ONLY: Laytypupw
*      USE GLOBAL,      ONLY:Iout,Hnew,Ibound
*      USE GWFBASMODULE, ONLY:HNOFLO
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION, EXTERNAL :: GW_func
*!     -----------------------------------------------------------------
*!     ARGUMENTS
*!     -----------------------------------------------------------------
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION rms, ferr
*      INTEGER ichld, irhld, ilhld, icfld, irfld, ilfld
*      INTEGER jj, ic, ir, il
*!     -----------------------------------------------------------------   
*      rms = 0.0D0
*      Fflux = 0.0D0
*      DO jj = 1, Numactive
*        il = Diag(jj, 1)
*        ir = Diag(jj, 2)
*        ic = Diag(jj, 3)
*        IF ( IBOUND(ic,ir,il).GT.0 ) THEN
*          IF( Laytypupw(il).GT.0 ) THEN
*            CALL TEMPFILLUN(ic, ir, il)
*          ELSE
*            CALL TEMPFILLCON(ic, ir, il)
*          END IF
*          ferr = GW_func(ic, ir, il)
*      
*          rms = rms + ferr**2.0D0
*          IF ( abs(ferr).GT.abs(Fflux) ) THEN
*            fflux = ferr
*            icfld = ic
*            irfld = ir
*            ilfld = il
*          END IF
*        END IF
*      END DO
*      RMS_func = rms**0.5D0
*      END FUNCTION RMS_func
!
!
!     -----------------------------------------------------------------
!     Calculates derivatives of conductance.
*      SUBROUTINE Dcon(kkiter)
*      USE GWFNWTMODULE
*      USE GWFUPWMODULE, ONLY:LAYTYPUPW
*      USE GLOBAL,      ONLY:Iout,HNEW,BOTM,LBOTM,NLAY,CC,CR,CV,NCOL,
*     +                      NROW
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION, EXTERNAL :: Dvert, Dhoriz
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION hh, ttop, bbot, Dv, Dh
*      INTEGER ij, ic, ir, il, i, kkiter
*!     -----------------------------------------------------------------   
*      DO ij = 1, Numactive
*        il = Diag(ij, 1)
*        ir = Diag(ij, 2)
*        ic = Diag(ij, 3)
*        DO I = 1,6
*          Dc(ij,I) = 0.0D0
*        END DO
*        hh = HNEW(ic,ir,il)
*        ttop = dble(BOTM(ic,ir,LBOTM(il)-1))
*        bbot = dble(BOTM(ic,ir,LBOTM(il)))
*        IF ( LAYTYPUPW(il)==0 ) THEN
*          Dv = 0.0D0
*          Dh = 0.0D0
*        ELSE
*          Dv = Dvert(hh,ttop,bbot)
*          Dh = (ttop-bbot)*Dhoriz(hh,ttop,bbot,il)
*          IF ( il.GT.1 ) THEN
*            IF ( hh.GT.HNEW(ic,ir,il-1)) Dc(ij,1) = 
*     +                                     Cv(ic,ir,il-1)*Dv
*          END IF
*          IF ( IR.GT.1 ) THEN
*            IF ( hh.GT.HNEW(ic,ir-1,il)) Dc(ij,2) = 
*     +                                     Cc(ic,ir-1,il)*Dh
*          END IF
*          IF ( IC.GT.1 ) THEN
*            IF ( hh.GT.HNEW(ic-1,ir,il)) Dc(ij,3) = 
*     +                                     Cr(ic-1,ir,il)*Dh
*          END IF
*          IF ( IC.LT.NCOL ) THEN
*            IF ( hh.GT.HNEW(ic+1,ir,il)) Dc(ij,4) = Cr(ic,ir,il)*Dh
*          END IF
*          IF ( IR.LT.NROW ) THEN
*            IF ( hh.GT.HNEW(ic,ir+1,il)) Dc(ij,5) = Cc(ic,ir,il)*Dh
*          END IF
*          IF ( IL.LT.NLAY ) THEN
*            IF ( hh.GT.HNEW(ic,ir,il+1)) Dc(ij,6) = Cv(ic,ir,il)*Dv
*          END IF
*        END IF
*      END DO
*      RETURN
*      END SUBROUTINE
!
!     -----------------------------------------------------------------
!     Returns sum of saturated thicknesses for cells connected
!     to ic,ir,il.
!     and head residuals
*      DOUBLE PRECISION FUNCTION Sum_sat(sum,ic,ir,il)
*      USE GWFNWTMODULE, ONLY:Icell
*      USE GWFUPWMODULE, ONLY:Sn
*      USE GLOBAL,       ONLY:Hnew,Ibound,Botm,Lbotm,Ncol,Nrow,Nlay
*      USE GWFBASMODULE, ONLY:HNOFLO
*      IMPLICIT NONE
*!     -----------------------------------------------------------------
*!     ARGUMENTS
*      INTEGER ic,ir,il,jj
*!     -----------------------------------------------------------------
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION sum, sat
*
*!     -----------------------------------------------------------------  
*      Sum_sat = 0.0D0
*      sat = 0.0D0
*      IF ( ic.GT.Ncol .OR. ic.LT.1 ) RETURN
*      IF ( ir.GT.Nrow .OR. ir.LT.1 ) RETURN
*      IF ( il.GT.Nlay .OR. il.LT.1 ) RETURN
*      jj = Icell(ic, ir, il)
*      IF ( jj.GT.0 ) sat = Sn(jj)
*      Sum_sat = sum + sat     
*      END FUNCTION Sum_sat
!
*      SUBROUTINE Jacobian(kkiter,kper,kstp)
*      USE GLOBAL, ONLY:Ncol, Nrow, Nlay, Ibound, Iout, Hnew, Botm,
*     +                 Lbotm, CR, CC, CV
*      USE GWFNWTMODULE
*      USE GWFUPWMODULE, ONLY: LAYTYPUPW
*      IMPLICIT NONE
*!     ------------------------------------------------------------------
*!     SPECIFICATIONS:
*!     ------------------------------------------------------------------
*      DOUBLE PRECISION, EXTERNAL :: GW_func, DVERT, DHORIZ
*      EXTERNAL TEMPFILLUN,TEMPFILLCON
*      INTEGER kkiter, kper
*!     -----------------------------------------------------------------
*!     LOCAL VARIABLES
*!     -----------------------------------------------------------------
*      DOUBLE PRECISION dum, ferr, Adiag, DD, tolf2
*      DOUBLE PRECISION zero, ttop, bbot, Dv, Dh, sum, tolf, botcheck
*      DOUBLE PRECISION term1, term2, term3, term4, term5, term6, coef
*      INTEGER ic, ir, il, icc, irr, ill
*      INTEGER ij, iltyp, kstp, I, I1, I2
*!     -----------------------------------------------------------------
*      zero = 0.0D0
*! Set conductance derivatives.
*      CALL Dcon(kkiter)
*! Fill Jacobian
*      DO ij = 1, Numactive
*        il = Diag(ij, 1)
*        ir = Diag(ij, 2)
*        ic = Diag(ij, 3)
*        Hchange(ij) = HNEW(ic,ir,il)
*        botcheck = dble(BOTM(ic,ir,il))
*! Constant head cells.
*        IF ( IBOUND(ic,ir,il).LT.0 ) THEN
*          A(IA(ij)) = 1.0D0
*          BB(ij) = HNEW(ic,ir,il)
*          DO I = IA(ij)+1,IA(ij+1)-1
*            A(I) = 0.0D0
*          END DO
*! Variable head cells.
*        ELSE
*          iltyp = LAYTYPUPW(il)
*          dum = 1.0D0
*! DIAGONAL FIRST
*          IF ( iltyp.GT.0 ) THEN
*            CALL TEMPFILLUN(ic, ir, il)
*          ELSE 
*            CALL TEMPFILLCON(ic, ir, il)
*          END IF
*          IF ( iltyp.GT.0 ) THEN
*            DO I = IA(ij)+1,IA(ij+1)-1
*              ill = Diag(JA(I), 1)
*              irr = Diag(JA(I), 2)
*              icc = Diag(JA(I), 3)
*              dd = 0.0D0
*              IF ( ill.LT.il ) THEN
*                dd = Dc(ij,1)
*              ELSEIF ( irr.LT.ir ) THEN
*                dd = Dc(ij,2)
*              ELSEIF ( icc.LT.ic ) THEN
*                dd = Dc(ij,3)
*              ELSEIF ( icc.GT.ic ) THEN
*                dd = Dc(ij,4)
*              ELSEIF ( irr.GT.ir ) THEN
*                dd = Dc(ij,5)
*              ELSEIF ( ill.GT.il ) THEN
*                dd = Dc(ij,6)
*              END IF
*                IF (abs(dd).GT.CLOSEZERO) A(IA(ij)) = A(IA(ij)) + dd*
*     +                         (HNEW(icc,irr,ill)-HNEW(ic,ir,il))
*            END DO
*          END IF
*          I = IA(ij)
*          A(I) = A(I)-Cvm1-Ccm1-Crm1-Cvv-Ccc-Crr+Hcoff
*! Calculate the right hand side
*          ferr = -GW_func(ic, ir, il)
*          BB(ij) = ferr
*          BB(ij) = BB(ij) + A(I)*HNEW(ic,ir,il)
*          sum  = 0.0D0
*          DO I = IA(ij)+1,IA(ij+1)-1
*            ill = Diag(JA(I), 1)
*            irr = Diag(JA(I), 2)
*            icc = Diag(JA(I), 3)
*            IF ( ill.LT.il ) THEN
*              A(I) = (Cvm1 + Dc(Icell(ic,ir,il-1),6)*(Hvm1-H))
*              BB(ij) = BB(ij) + A(I)*Hvm1
*              sum = sum + A(I)
*            ELSEIF ( irr.LT.ir ) THEN
*              A(I) = (Ccm1 + Dc(Icell(ic,ir-1,il),5)*(Hrm1-H))
*              BB(ij) = BB(ij) + A(I)*Hrm1
*              sum = sum + A(I)
*            ELSEIF ( icc.LT.ic ) THEN
*              A(I) = (Crm1 + Dc(Icell(ic-1,ir,il),4)*(Hcm1-H))
*              BB(ij) = BB(ij) + A(I)*Hcm1
*              sum = sum + A(I)
*            ELSEIF ( icc.GT.ic ) THEN
*              A(I) = (Crr + Dc(Icell(ic+1,ir,il),3)*(Hcp1-H))
*              BB(ij) = BB(ij) + A(I)*Hcp1
*              sum = sum + A(I)
*            ELSEIF ( irr.GT.ir ) THEN
*              A(I) = (Ccc + Dc(Icell(ic,ir+1,il),2)*(Hrp1-H))
*              BB(ij) = BB(ij) + A(I)*Hrp1
*              sum = sum + A(I)
*            ELSEIF ( ill.GT.il ) THEN
*              A(I) = (Cvv + Dc(Icell(ic,ir,il+1),1)*(Hvp1-H))
*              BB(ij) = BB(ij) + A(I)*Hvp1
*              sum = sum + A(I)
*            ENDIF
*          END DO
*! Check for row values bad for linear solver
*          Adiag = ABS(A(IA(ij)))
*          tolf = 1.0D-9
*          IF ( Adiag.LT.tolf ) THEN
*            A(IA(ij)) = 1.0D0
*           BB(ij) = BB(ij) + HNEW(ic,ir,il)
*          ELSE IF ( abs(sum).LT.tolf ) THEN
*            A(IA(ij)) = A(IA(ij)) + tolf
*            BB(ij) = BB(ij) + tolf*HNEW(ic,ir,il)
*          END IF
*        END IF
*!        I1 = IA(ij)
*!        I2 = IA(ij+1)-1
*!        do i=i1,i2
*!        if ( kstp.eq.2)then
*!        if(ic==166.and.ir==159.and.il==1)then
*!        WRITE(IOUT,66)ij,BB(ij),HNEW(ic,ir,il),BOTM(ic,ir,il-1),
*!     +                BOTM(ic,ir,il),HCOFF,RHSS,(A(I),I=I1,I2)
*!       write(iout,*)IA(ij),i,A(I),BB(ij)
*!        end if
*!!        end if      
*! 66      FORMAT(I9,1X,3G15.6,2X,11G15.6)
*!       end do
*      END DO
*      RETURN
*      END SUBROUTINE Jacobian
*!
*      
*
      end module unusedcode
