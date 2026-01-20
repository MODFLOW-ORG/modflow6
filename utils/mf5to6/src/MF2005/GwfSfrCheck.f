C
C
C     ******************************************************************
C     CHECK FOR STEAMBED BELOW CELL BOTTOM. RECORD REACHES FOR PRINTING
C     ******************************************************************
      MODULE GwfSfrCheckModule

      USE GWFSFRMODULE,ONLY:ISTRM,STRM,NSTRM
      USE GLOBAL,ONLY:BOTM,IBOUND,LBOTM
      use SimPHMFModule, only: ustop
      implicit none
      type check_bot
        integer ltype,irchnum,iflag,iunit
      end type check_bot
      public check_bot

      CONTAINS

      FUNCTION ICHKSTRBOT(self)
      type (check_bot), intent(in) :: self
      INTEGER JRCH,IRCH,KRCH,JSEG,ISEG,ICHKSTRBOT
      ICHKSTRBOT = 0
      KRCH = ISTRM(1,self%IRCHNUM)
      IRCH = ISTRM(2,self%IRCHNUM)
      JRCH = ISTRM(3,self%IRCHNUM)
      JSEG = ISTRM(4,self%IRCHNUM)
      ISEG = ISTRM(5,self%IRCHNUM)
      IF ( self%LTYPE.GT.0  .AND. IBOUND(JRCH,IRCH,KRCH).GT.0 ) THEN
        IF ( STRM(4, self%IRCHNUM)-BOTM(JRCH,IRCH,LBOTM(KRCH))
     +                                      .LT.-1.0E-12 ) THEN
          IF ( self%IFLAG.EQ.0 ) THEN
          WRITE(self%IUNIT,*)
          WRITE(self%IUNIT,*)' REACHES WITH ALTITUDE ERRORS:'
          WRITE(self%IUNIT,*)'   LAY    ROW    COL    SEG  REACH      ',
     +                'STR.ELEV.      CELL-BOT.'
          END IF
          WRITE(self%IUNIT,100)KRCH,IRCH,JRCH,JSEG,ISEG,
     +                STRM(4, self%IRCHNUM),BOTM(JRCH,IRCH,LBOTM(KRCH))
          ICHKSTRBOT = 1
        END IF
      END IF
      IF ( self%IFLAG.GT.0 .AND. self%IRCHNUM.EQ.NSTRM ) THEN
        WRITE(self%IUNIT,*)' MODEL STOPPING DUE TO REACH ALTITUDE ERROR'
        CALL USTOP(' ')
      END IF
  100 FORMAT(5I7,2F15.7)
      END FUNCTION ICHKSTRBOT

      END MODULE GwfSfrCheckModule
C
