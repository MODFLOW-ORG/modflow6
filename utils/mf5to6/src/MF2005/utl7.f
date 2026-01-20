! ERB 1/20/2015 -- Convert to module to make explicit interfaces
!                  and allow optional arguments.
      module utl7module

        use ConstantsModule, only: LENOBSNAME, LINELENGTH
        use ConstantsPHMFModule, only: LENOBSNAMENEW
        use FileTypeModule, only: FileType
        use FileListModule, only: FileListType
        use OpenSpecModule, only: ACCESS, ACTION, FORM
        use PRECUTLSMOD, only: GET_BINARY_HEAD_DATASET, HEADPRECISION, 
     &                         GET_BINARY_FILE_LENGTH
        use SimPHMFModule, only: store_error, ustop
        
        private
        public :: URWORD, UPCASE, U1DREL, U2DINT, U2DREL, USTOPx
!        public :: UBDSV2, UBDSV3, UBDSV4, UBDSVA, u
        public :: ULAPRS, ULAPRW, ULAPRWC, ULSTLB, ULSTRD, UMESPR
        public :: URDCOM, assign_ncharsizes_flow, build_obsname
        public :: UBUDSV, UBDSVB, UBDSV3

      contains

      SUBROUTINE URWORD(LINE,ICOL,ISTART,ISTOP,NCODE,N,R,IOUT,IN)
C     ******************************************************************
C     ROUTINE TO EXTRACT A WORD FROM A LINE OF TEXT, AND OPTIONALLY
C     CONVERT THE WORD TO A NUMBER.
C        ISTART AND ISTOP WILL BE RETURNED WITH THE STARTING AND
C          ENDING CHARACTER POSITIONS OF THE WORD.
C        THE LAST CHARACTER IN THE LINE IS SET TO BLANK SO THAT IF ANY
C          PROBLEMS OCCUR WITH FINDING A WORD, ISTART AND ISTOP WILL
C          POINT TO THIS BLANK CHARACTER.  THUS, A WORD WILL ALWAYS BE
C          RETURNED UNLESS THERE IS A NUMERIC CONVERSION ERROR.  BE SURE
C          THAT THE LAST CHARACTER IN LINE IS NOT AN IMPORTANT CHARACTER
C          BECAUSE IT WILL ALWAYS BE SET TO BLANK.
C        A WORD STARTS WITH THE FIRST CHARACTER THAT IS NOT A SPACE OR
C          COMMA, AND ENDS WHEN A SUBSEQUENT CHARACTER THAT IS A SPACE
C          OR COMMA.  NOTE THAT THESE PARSING RULES DO NOT TREAT TWO
C          COMMAS SEPARATED BY ONE OR MORE SPACES AS A NULL WORD.
C        FOR A WORD THAT BEGINS WITH "'", THE WORD STARTS WITH THE
C          CHARACTER AFTER THE QUOTE AND ENDS WITH THE CHARACTER
C          PRECEDING A SUBSEQUENT QUOTE.  THUS, A QUOTED WORD CAN
C          INCLUDE SPACES AND COMMAS.  THE QUOTED WORD CANNOT CONTAIN
C          A QUOTE CHARACTER.
C        IF NCODE IS 1, THE WORD IS CONVERTED TO UPPER CASE.
C        IF NCODE IS 2, THE WORD IS CONVERTED TO AN INTEGER.
C        IF NCODE IS 3, THE WORD IS CONVERTED TO A REAL NUMBER.
C        NUMBER CONVERSION ERROR IS WRITTEN TO UNIT IOUT IF IOUT IS
C          POSITIVE; ERROR IS WRITTEN TO DEFAULT OUTPUT IF IOUT IS 0;
C          NO ERROR MESSAGE IS WRITTEN IF IOUT IS NEGATIVE.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      double precision :: r
      CHARACTER*(*) LINE
      CHARACTER*20 STRING
      CHARACTER*30 RW
      CHARACTER*1 TAB
      character(len=LINELENGTH) :: msg
C     ------------------------------------------------------------------
      TAB=CHAR(9)
C
C1------Set last char in LINE to blank and set ISTART and ISTOP to point
C1------to this blank as a default situation when no word is found.  If
C1------starting location in LINE is out of bounds, do not look for a
C1------word.
      LINLEN=LEN(LINE)
      LINE(LINLEN:LINLEN)=' '
      ISTART=LINLEN
      ISTOP=LINLEN
      LINLEN=LINLEN-1
      IF(ICOL.LT.1 .OR. ICOL.GT.LINLEN) GO TO 100
C
C2------Find start of word, which is indicated by first character that
C2------is not a blank, a comma, or a tab.
      DO 10 I=ICOL,LINLEN
      IF(LINE(I:I).NE.' ' .AND. LINE(I:I).NE.','
     &    .AND. LINE(I:I).NE.TAB) GO TO 20
10    CONTINUE
      ICOL=LINLEN+1
      GO TO 100
C
C3------Found start of word.  Look for end.
C3A-----When word is quoted, only a quote can terminate it.
20    IF(LINE(I:I).EQ.'''') THEN
         I=I+1
         IF(I.LE.LINLEN) THEN
            DO 25 J=I,LINLEN
            IF(LINE(J:J).EQ.'''') GO TO 40
25          CONTINUE
         END IF
C
C3B-----When word is not quoted, space, comma, or tab will terminate.
      ELSE
         DO 30 J=I,LINLEN
         IF(LINE(J:J).EQ.' ' .OR. LINE(J:J).EQ.','
     &    .OR. LINE(J:J).EQ.TAB) GO TO 40
30       CONTINUE
      END IF
C
C3C-----End of line without finding end of word; set end of word to
C3C-----end of line.
      J=LINLEN+1
C
C4------Found end of word; set J to point to last character in WORD and
C-------set ICOL to point to location for scanning for another word.
40    ICOL=J+1
      J=J-1
      IF(J.LT.I) GO TO 100
      ISTART=I
      ISTOP=J
C
C5------Convert word to upper case and RETURN if NCODE is 1.
      IF(NCODE.EQ.1) THEN
         IDIFF=ICHAR('a')-ICHAR('A')
         DO 50 K=ISTART,ISTOP
            IF(LINE(K:K).GE.'a' .AND. LINE(K:K).LE.'z')
     1             LINE(K:K)=CHAR(ICHAR(LINE(K:K))-IDIFF)
50       CONTINUE
         RETURN
      END IF
C
C6------Convert word to a number if requested.
100   IF(NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
         RW=' '
         L=30-ISTOP+ISTART
         IF(L.LT.1) GO TO 200
         RW(L:30)=LINE(ISTART:ISTOP)
         IF(NCODE.EQ.2) READ(RW,'(I30)',ERR=200) N
         IF(NCODE.EQ.3) READ(RW,'(F30.0)',ERR=200) R
      END IF
      RETURN
C
C7------Number conversion error.
200   IF(NCODE.EQ.3) THEN
         STRING= 'A REAL NUMBER'
         L=13
      ELSE
         STRING= 'AN INTEGER'
         L=10
      END IF
C
C7A-----If output unit is negative, set last character of string to 'E'.
      IF(IOUT.LT.0) THEN
         N=0
         R=0.
         LINE(LINLEN+1:LINLEN+1)='E'
         RETURN
C
C7B-----If output unit is positive; write a message to output unit.
      ELSE IF(IOUT.GT.0) THEN
         IF(IN.GT.0) THEN
            WRITE(IOUT,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
            write(msg,301)IN,LINE(ISTART:ISTOP),STRING(1:L)
            call store_error(msg)
            write(msg,311)trim(line)
            call store_error(msg)
         ELSE
            WRITE(IOUT,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
            write(msg,302)LINE(ISTART:ISTOP),STRING(1:L)
            call store_error(msg)
            write(msg,311)trim(line)
            call store_error(msg)
         END IF
201      FORMAT(1X,/1X,'FILE UNIT ',I4,' : ERROR CONVERTING "',A,
     1       '" TO ',A,' IN LINE:',/1X,A)
202      FORMAT(1X,/1X,'KEYBOARD INPUT : ERROR CONVERTING "',A,
     1       '" TO ',A,' IN LINE:',/1X,A)
301      FORMAT('FILE UNIT ',I4,' : ERROR CONVERTING "',A,
     &       '" TO ',A,' IN LINE:')
311      format(a)
302      FORMAT('KEYBOARD INPUT : ERROR CONVERTING "',A,
     &       '" TO ',A,' IN LINE:')
C
C7C-----If output unit is 0; write a message to default output.
      ELSE
         IF(IN.GT.0) THEN
            WRITE(*,201) IN,LINE(ISTART:ISTOP),STRING(1:L),LINE
         ELSE
            WRITE(*,202) LINE(ISTART:ISTOP),STRING(1:L),LINE
         END IF
      END IF
C
C7D-----STOP after writing message.
      CALL USTOP(' ')
      END SUBROUTINE URWORD
C*******************************************************************************
      SUBROUTINE UPCASE(WORD)
C     ******************************************************************
C     CONVERT A CHARACTER STRING TO ALL UPPER CASE
C     ******************************************************************
C       SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER WORD*(*)
C
C1------Compute the difference between lowercase and uppercase.
      L = LEN(WORD)
      IDIFF=ICHAR('a')-ICHAR('A')
C
C2------Loop through the string and convert any lowercase characters.
      DO 10 K=1,L
      IF(WORD(K:K).GE.'a' .AND. WORD(K:K).LE.'z')
     1   WORD(K:K)=CHAR(ICHAR(WORD(K:K))-IDIFF)
10    CONTINUE
C
C3------return.
      RETURN
      END SUBROUTINE UPCASE
C*******************************************************************************
      SUBROUTINE URDCOM(IN,IOUT,LINE)
C     ******************************************************************
C     READ COMMENTS FROM A FILE AND PRINT THEM.  RETURN THE FIRST LINE
C     THAT IS NOT A COMMENT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*(*) LINE
C     ------------------------------------------------------------------
C
C1------Read a line
   10 READ(IN,'(A)') LINE
C
C2------If the line does not start with "#", return.
      IF(LINE(1:1).NE.'#') RETURN
C
C3------Find the last non-blank character.
      L=LEN(LINE)
      DO 20 I=L,1,-1
      IF(LINE(I:I).NE.' ') GO TO 30
   20 CONTINUE
C
C4------Print the line up to the last non-blank character if IOUT>0.
   30 IF (IOUT.GT.0) WRITE(IOUT,'(1X,A)') LINE(1:I)
      GO TO 10
C
      END SUBROUTINE URDCOM
C*******************************************************************************
      SUBROUTINE ULSTRD(NLIST,RLIST,LSTBEG,LDIM,MXLIST,IAL,INPACK,IOUT,
     1     LABEL,CAUX,NCAUX,NAUX,IFREFM,NCOL,NROW,NLAY,ISCLOC1,ISCLOC2,
     2     IPRFLG)
C     ******************************************************************
C     Read and print a list.  NAUX of the values in the list are
C     optional -- auxiliary data.
C     ******************************************************************
      CHARACTER*(*) LABEL
      CHARACTER*16 CAUX(NCAUX)
      double precision :: RLIST(LDIM,MXLIST), r, sfac
      CHARACTER*200 LINE,FNAME
      CHARACTER*20 FMTARG, ACCARG
      CHARACTER*20 FILEFMT
      CHARACTER*30 CERR
      LOGICAL LVAL
      DATA NUNOPN/99/
C     ------------------------------------------------------------------
C
C1------If the list is empty, return.
      IF (NLIST.EQ.0) RETURN
C
C2------Check for and decode EXTERNAL and OPEN/CLOSE records.
      IN=INPACK
      ICLOSE=0
      IBINARY=0
      READ(IN,'(A)') LINE
      SFAC=1.
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
      IF(LINE(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
         IN=I
C          TEST IF EXTERNAL FILE IS OPEN
         INQUIRE( UNIT=IN, OPENED=LVAL )
         IF ( LVAL.EQV. .FALSE. ) THEN
           WRITE ( CERR,110 ) IN
  110    FORMAT('External unit ', I4,' is not open')
           WRITE ( IOUT,'(1X,A)' ) CERR
           CALL USTOPx(CERR)
         END IF
C          TEST IF OPEN EXTERNAL FILE IS FORMATTED OR UNFORMATTED/BINARY
C          SEE FORM VARIABLE IN OpenSpecModule
         FMTARG=FORM
         INQUIRE( UNIT=IN, FORM=FILEFMT )
         IF ( FILEFMT.EQ.FMTARG ) THEN
           IBINARY=1
         END IF
         IF(IPRFLG.EQ.1) THEN
           IF ( IBINARY.NE.1 ) THEN
             WRITE(IOUT,111) IN
  111        FORMAT(1X,'Reading list on unit ',I4)
           ELSE
             WRITE(IOUT,1111) IN
 1111        FORMAT(1X,'Reading list on binary unit ',I4)
           END IF
         END IF
         IF (IBINARY.NE.1) READ(IN,'(A)') LINE
      ELSE IF(LINE(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=LINE(ISTART:ISTOP)
         IN=NUNOPN
C          TEST IF OPEN\CLOSE FILE EXISTS
         INQUIRE( FILE=FNAME, EXIST=LVAL )
         IF ( LVAL.EQV. .FALSE. ) THEN
           WRITE ( IOUT,112 ) LINE(ISTART:ISTOP)
  112      FORMAT('Specified OPEN/CLOSE file ',(A),' does not exit')
           CALL USTOPx('Specified OPEN/CLOSE file does not exit')
         END IF
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,0,N,R,IOUT,IN)
         IF (LINE(ISTART:ISTOP).EQ.'(BINARY)') THEN
           IBINARY=1
           IF(IPRFLG.EQ.1)WRITE(IOUT,1115) IN,FNAME
 1115      FORMAT(1X,/1X,'OPENING BINARY FILE ON UNIT ',I4,':',/1X,A)
           FMTARG=FORM
           ACCARG=ACCESS
           OPEN(UNIT=IN,FILE=FNAME,ACTION=ACTION(1),FORM=FMTARG,
     1          ACCESS=ACCARG,STATUS='OLD')
         ELSE
           IF(IPRFLG.EQ.1)WRITE(IOUT,115) IN,FNAME
  115      FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
           OPEN(UNIT=IN,FILE=FNAME,ACTION=ACTION(1))
         END IF
         ICLOSE=1
         IF (IBINARY.NE.1) READ(IN,'(A)') LINE
      END IF
C
C3------Check for SFAC record in ascii records.
      IF (IBINARY.NE.1) THEN
        LLOC=1
        CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,I,R,IOUT,IN)
        IF(LINE(ISTART:ISTOP).EQ.'SFAC') THEN
           CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,I,SFAC,IOUT,IN)
           IF(IPRFLG.EQ.1) THEN
             WRITE(IOUT,116) SFAC
  116        FORMAT(1X,'LIST SCALING FACTOR=',1PG12.5)
             IF(ISCLOC1.EQ.ISCLOC2) THEN
                WRITE(IOUT,113) ISCLOC1
  113         FORMAT(1X,'(THE SCALE FACTOR WAS APPLIED TO FIELD',I2,')')
             ELSE
                WRITE(IOUT,114) ISCLOC1,ISCLOC2
  114           FORMAT(1X,'(THE SCALE FACTOR WAS APPLIED TO FIELDS',
     1             I2,'-',I2,')')
             END IF
           ENDIF
           READ(IN,'(A)') LINE
        END IF
      END IF
C
C3------Write a label for the list if the list will be printed.
      IF(IPRFLG.EQ.1) THEN
         WRITE(IOUT,'(1X)')
         CALL ULSTLB(IOUT,LABEL,CAUX,NCAUX,NAUX)
      END IF
C
C4------Setup indices for reading the list
      NREAD2=LDIM-IAL
      NREAD1=NREAD2-NAUX
      N=NLIST+LSTBEG-1
C
C4A-----READ THE LIST -- BINARY OR ASCII
      IF (IBINARY.NE.0) THEN
        READ(IN) ((RLIST(JJ,II),JJ=1,NREAD2),II=LSTBEG,N)
      ELSE
C
C5------READ AN ASCII LIST
        DO 240 II=LSTBEG,N
C
C5A-----Read a line into the buffer.  (The first line has already been
C5A-----read to scan for EXTERNAL and SFAC records.)
        IF(II.NE.LSTBEG) READ(IN,'(A)') LINE
C
C5B-----Get the non-optional values from the line.
        IF(IFREFM.EQ.0) THEN
          READ(LINE,'(3I10,9F10.0)') K,I,J,(RLIST(JJ,II),JJ=4,NREAD1)
           LLOC=10*NREAD1+1
        ELSE
          LLOC=1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,K,R,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,I,R,IOUT,IN)
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,J,R,IOUT,IN)
          DO 200 JJ=4,NREAD1
          CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,RLIST(JJ,II),
     1                IOUT,IN)
200       CONTINUE
        END IF
        RLIST(1,II)=K
        RLIST(2,II)=I
        RLIST(3,II)=J
C
C5E-----Get the optional values from the line
        IF(NAUX.GT.0) THEN
           DO 210 JJ=NREAD1+1,NREAD2
           CALL URWORD(LINE,LLOC,ISTART,ISTOP,3,IDUM,RLIST(JJ,II),
     1                 IOUT,IN)
210        CONTINUE
        END IF
C
240     CONTINUE
      ENDIF
C
C6----SCALE THE DATA AND CHECK
      DO 250 II=LSTBEG,N
      K=RLIST(1,II)
      I=RLIST(2,II)
      J=RLIST(3,II)
C
C6A------Scale fields ISCLOC1-ISCLOC2 by SFAC
      DO 204 ILOC=ISCLOC1,ISCLOC2
        RLIST(ILOC,II)=RLIST(ILOC,II)*SFAC
204   CONTINUE
C
C6B-----Write the values that were read if IPRFLG is 1.
      NN=II-LSTBEG+1
      IF(IPRFLG.EQ.1)
     1    WRITE(IOUT,205) NN,K,I,J,(RLIST(JJ,II),JJ=4,NREAD2)
205       FORMAT(1X,I6,I7,I7,I7,26G16.4)
C
C6C-----Check for illegal grid location
      IF(K.LT.1 .OR. K.GT.NLAY) THEN
         WRITE(IOUT,*) ' Layer number in list is outside of the grid'
         CALL USTOPx(' ')
      END IF
      IF(I.LT.1 .OR. I.GT.NROW) THEN
         WRITE(IOUT,*) ' Row number in list is outside of the grid'
         CALL USTOPx(' ')
      END IF
      IF(J.LT.1 .OR. J.GT.NCOL) THEN
         WRITE(IOUT,*) ' Column number in list is outside of the grid'
         CALL USTOPx(' ')
      END IF
  250 CONTINUE
C
C7------Done reading the list.  If file is open/close, close it.
      IF(ICLOSE.NE.0) CLOSE(UNIT=IN)
C
      RETURN
      END SUBROUTINE ULSTRD
C*******************************************************************************
      SUBROUTINE ULSTLB(IOUT,LABEL,CAUX,NCAUX,NAUX)
C     ******************************************************************
C     PRINT A LABEL FOR A LIST
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*(*) LABEL
      CHARACTER*16 CAUX(NCAUX)
      CHARACTER*400 BUF
      CHARACTER*1 DASH(400)
      DATA DASH/400*'-'/
C     ------------------------------------------------------------------
C
C1------Construct the complete label in BUF.  Start with BUF=LABEL.
      BUF=LABEL
C
C2------Add auxiliary data names if there are any.
      NBUF=LEN(LABEL)+9
      IF(NAUX.GT.0) THEN
         DO 10 I=1,NAUX
         N1=NBUF+1
         NBUF=NBUF+16
         BUF(N1:NBUF)=CAUX(I)
10       CONTINUE
      END IF
C
C3------Write the label.
      WRITE(IOUT,103) BUF(1:NBUF)
  103 FORMAT(1X,A)
C
C4------Add a line of dashes.
      WRITE(IOUT,104) (DASH(J),J=1,NBUF)
  104 FORMAT(1X,400A)
C
C5------Return.
      RETURN
      END SUBROUTINE ULSTLB
C*******************************************************************************
      SUBROUTINE U1DREL(A,ANAME,JJ,IN,IOUT,constantval,constant,iprnarg)
C     ******************************************************************
C     ROUTINE TO INPUT 1-D REAL DATA MATRICES
C       A IS ARRAY TO INPUT
C       ANAME IS 24 CHARACTER DESCRIPTION OF A
C       JJ IS NO. OF ELEMENTS
C       IN IS INPUT UNIT
C       IOUT IS OUTPUT UNIT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*24 ANAME
      double precision :: A(JJ), r, cnstnt
      double precision, optional :: constantval
      logical, optional :: constant
      integer, optional :: iprnarg
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*200 FNAME
      DATA NUNOPN/99/
C     ------------------------------------------------------------------
C
C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)') CNTRL
C
C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      if (present(constant)) constant = .false.
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
         if (present(constant)) constant = .true.
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=IN
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN)
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
         WRITE(IOUT,15) LOCAT,FNAME
   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         OPEN(UNIT=LOCAT,FILE=FNAME,ACTION=ACTION(1))
         ICLOSE=1
      ELSE
C
C2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
C2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=500) LOCAT,CNSTNT,FMTIN,IPRN
         if (present(constant) .and. locat==0) constant = .true.
         if (present(iprnarg)) then
           iprnarg = iprn
         endif
    1    FORMAT(I10,F10.0,A20,I10)
      END IF
C
C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,3,N,CNSTNT,IOUT,IN)
         IF(LOCAT.GT.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
            FMTIN=CNTRL(ISTART:ISTOP)
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
            if (present(iprnarg)) then
              iprnarg = iprn
            endif
         END IF
      END IF
      !
      if (present(constantval)) then
        constantval = cnstnt
      endif
C
C4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT.GT.0) GO TO 90
C
C4A-----LOCAT <0 OR =0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
      DO 80 J=1,JJ
   80 A(J)=CNSTNT
      WRITE(IOUT,3) ANAME,CNSTNT
    3 FORMAT(1X,/1X,A,' =',1P,G14.6)
      RETURN
C
C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
   90 CONTINUE
      WRITE(IOUT,5) ANAME,LOCAT,FMTIN
    5 FORMAT(1X,///11X,A,/
     1       1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A20)
      IF(FMTIN.EQ.'(FREE)') THEN
      READ(LOCAT,*) (A(J),J=1,JJ)
      ELSE
         READ(LOCAT,FMTIN) (A(J),J=1,JJ)
      END IF
      IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
C
C5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
      ZERO=0.
      IF(CNSTNT.EQ.ZERO) GO TO 120
      DO 100 J=1,JJ
  100 A(J)=A(J)*CNSTNT
C
C6------IF PRINT CODE (IPRN) =0 OR >0 THEN PRINT ARRAY VALUES.
120   CONTINUE
      IF(IPRN.EQ.0) THEN
         WRITE(IOUT,1001) (A(J),J=1,JJ)
1001     FORMAT((1X,1PG12.5,9(1X,G12.5)))
      ELSE IF(IPRN.GT.0) THEN
         WRITE(IOUT,1002) (A(J),J=1,JJ)
1002     FORMAT((1X,1PG12.5,4(1X,G12.5)))
      END IF
C
C7------RETURN
      RETURN
C
C8------CONTROL RECORD ERROR.
500   WRITE(IOUT,502) ANAME
502   FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
      WRITE(IOUT,'(1X,A)') CNTRL
      CALL USTOPx(' ')
      END SUBROUTINE U1DREL
C*******************************************************************************
      SUBROUTINE U2DINT(IA,ANAME,II,JJ,K,IN,IOUT)
C     ******************************************************************
C     ROUTINE TO INPUT 2-D INTEGER DATA MATRICES
C       IA IS ARRAY TO INPUT
C       ANAME IS 24 CHARACTER DESCRIPTION OF IA
C       II IS NO. OF ROWS
C       JJ IS NO. OF COLS
C       K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --
C              IF K=0, NO LAYER IS PRINTED
C              IF K<0, CROSS SECTION IS PRINTED)
C       IN IS INPUT UNIT
C       IOUT IS OUTPUT UNIT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*24 ANAME
      DIMENSION IA(JJ,II)
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*200 FNAME
      double precision :: r
      DATA NUNOPN/99/
C     ------------------------------------------------------------------
C
C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)') CNTRL
C
C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=IN
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN)
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
         WRITE(IOUT,15) LOCAT,FNAME
   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         ICLOSE=1
      ELSE
C
C2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
C2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=600) LOCAT,ICONST,FMTIN,IPRN
    1    FORMAT(I10,I10,A20,I10)
      END IF
C
C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,ICONST,R,IOUT,IN)
         IF(LOCAT.NE.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
            FMTIN=CNTRL(ISTART:ISTOP)
            IF(ICLOSE.NE.0) THEN
               IF(FMTIN.EQ.'(BINARY)') THEN
                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM=FORM,ACCESS=ACCESS,
     &                 ACTION=ACTION(1))
               ELSE
                  OPEN(UNIT=LOCAT,FILE=FNAME,ACTION=ACTION(1))
               END IF
            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'(BINARY)') LOCAT=-LOCAT
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
         END IF
      END IF
C
C4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT.EQ.0) THEN
C
C4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO ICONST. RETURN.
        DO 80 I=1,II
        DO 80 J=1,JJ
   80   IA(J,I)=ICONST
        IF(K.GT.0) WRITE(IOUT,82) ANAME,ICONST,K
   82   FORMAT(1X,/1X,A,' =',I15,' FOR LAYER',I4)
        IF(K.LE.0) WRITE(IOUT,83) ANAME,ICONST
   83   FORMAT(1X,/1X,A,' =',I15)
        RETURN
      ELSE IF(LOCAT.GT.0) THEN
C
C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
        IF(K.GT.0) THEN
           WRITE(IOUT,94) ANAME,K,LOCAT,FMTIN
   94      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
        ELSE IF(K.EQ.0) THEN
           WRITE(IOUT,95) ANAME,LOCAT,FMTIN
   95      FORMAT(1X,///11X,A,/
     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
        ELSE
           WRITE(IOUT,96) ANAME,LOCAT,FMTIN
   96      FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
        END IF
        DO 100 I=1,II
        IF(FMTIN.EQ.'(FREE)') THEN
           READ(LOCAT,*) (IA(J,I),J=1,JJ)
        ELSE
           READ(LOCAT,FMTIN) (IA(J,I),J=1,JJ)
        END IF
  100   CONTINUE
      ELSE
C
C4C-----LOCAT<0; READ UNFORMATTED RECORD CONTAINING ARRAY VALUES.
        LOCAT=-LOCAT
        IF(K.GT.0) THEN
           WRITE(IOUT,201) ANAME,K,LOCAT
  201      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
     1      1X,'READING BINARY ON UNIT ',I4)
        ELSE IF(K.EQ.0) THEN
           WRITE(IOUT,202) ANAME,LOCAT
  202      FORMAT(1X,///11X,A,/
     1      1X,'READING BINARY ON UNIT ',I4)
        ELSE
           WRITE(IOUT,203) ANAME,LOCAT
  203      FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
     1      1X,'READING BINARY ON UNIT ',I4)
        END IF
        READ(LOCAT)
        READ(LOCAT) IA
      END IF
C
C5------IF ICONST NOT ZERO THEN MULTIPLY ARRAY VALUES BY ICONST.
      IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
      IF(ICONST.EQ.0) GO TO 320
      DO 310 I=1,II
      DO 310 J=1,JJ
      IA(J,I)=IA(J,I)*ICONST
  310 CONTINUE
C
C6------IF PRINT CODE (IPRN) <0 THEN RETURN.
  320 IF(IPRN.LT.0) RETURN
C
C7------PRINT COLUMN NUMBERS AT TOP OF PAGE.
      IF(IPRN.GT.9 .OR. IPRN.EQ.0) IPRN=6
      GO TO(401,402,403,404,405,406,407,408,409), IPRN
401   CALL UCOLNO(1,JJ,4,60,2,IOUT)
      GO TO 500
402   CALL UCOLNO(1,JJ,4,40,3,IOUT)
      GO TO 500
403   CALL UCOLNO(1,JJ,4,30,4,IOUT)
      GO TO 500
404   CALL UCOLNO(1,JJ,4,25,5,IOUT)
      GO TO 500
405   CALL UCOLNO(1,JJ,4,20,6,IOUT)
      GO TO 500
406   CALL UCOLNO(1,JJ,4,10,12,IOUT)
      GO TO 500
407   CALL UCOLNO(1,JJ,4,25,3,IOUT)
      GO TO 500
408   CALL UCOLNO(1,JJ,4,15,5,IOUT)
      GO TO 500
409   CALL UCOLNO(1,JJ,4,10,7,IOUT)
C
C8------PRINT EACH ROW IN THE ARRAY.
500   DO 510 I=1,II
      GO TO(501,502,503,504,505,506,507,508,509), IPRN
C
C----------------FORMAT 60I1
  501 WRITE(IOUT,551) I,(IA(J,I),J=1,JJ)
  551 FORMAT(1X,I3,1X,60(1X,I1):/(5X,60(1X,I1)))
      GO TO 510
C
C----------------FORMAT 40I2
  502 WRITE(IOUT,552) I,(IA(J,I),J=1,JJ)
  552 FORMAT(1X,I3,1X,40(1X,I2):/(5X,40(1X,I2)))
      GO TO 510
C
C----------------FORMAT 30I3
  503 WRITE(IOUT,553) I,(IA(J,I),J=1,JJ)
  553 FORMAT(1X,I3,1X,30(1X,I3):/(5X,30(1X,I3)))
      GO TO 510
C
C----------------FORMAT 25I4
  504 WRITE(IOUT,554) I,(IA(J,I),J=1,JJ)
  554 FORMAT(1X,I3,1X,25(1X,I4):/(5X,25(1X,I4)))
      GO TO 510
C
C----------------FORMAT 20I5
  505 WRITE(IOUT,555) I,(IA(J,I),J=1,JJ)
  555 FORMAT(1X,I3,1X,20(1X,I5):/(5X,20(1X,I5)))
      GO TO 510
C
C----------------FORMAT 10I11
  506 WRITE(IOUT,556) I,(IA(J,I),J=1,JJ)
  556 FORMAT(1X,I3,1X,10(1X,I11):/(5X,10(1X,I11)))
      GO TO 510
C
C----------------FORMAT 25I2
  507 WRITE(IOUT,557) I,(IA(J,I),J=1,JJ)
  557 FORMAT(1X,I3,1X,25(1X,I2):/(5X,25(1X,I2)))
      GO TO 510
C
C----------------FORMAT 15I4
  508 WRITE(IOUT,558) I,(IA(J,I),J=1,JJ)
  558 FORMAT(1X,I3,1X,15(1X,I4):/(5X,15(1X,I4)))
      GO TO 510
C
C----------------FORMAT 10I6
  509 WRITE(IOUT,559) I,(IA(J,I),J=1,JJ)
  559 FORMAT(1X,I3,1X,10(1X,I6):/(5X,10(1X,I6)))
C
  510 CONTINUE
C
C9------RETURN
      RETURN
C
C10-----CONTROL RECORD ERROR.
  600 IF(K.GT.0) THEN
         WRITE(IOUT,601) ANAME,K
  601    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,
     1     ' FOR LAYER',I4,':')
      ELSE
         WRITE(IOUT,602) ANAME
  602    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
      END IF
      WRITE(IOUT,'(1X,A)') CNTRL
      CALL USTOPx(' ')
      END SUBROUTINE U2DINT
C*******************************************************************************
      SUBROUTINE U2DREL(A,ANAME,II,JJ,K,IN,IOUT,constantval,constant,
     &    iprnarg,Mf2005Files)
      ! To support reading of input from an external, non-OPEN/CLOSE file,
      ! the Mf2005Files argument must be present and pointing (ultimately)
      ! to model%Mf2005Files.
C     ******************************************************************
C     ROUTINE TO INPUT 2-D REAL DATA MATRICES
C       A IS ARRAY TO INPUT
C       ANAME IS 24 CHARACTER DESCRIPTION OF A
C       II IS NO. OF ROWS
C       JJ IS NO. OF COLS
C       K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --)
C              IF K=0, NO LAYER IS PRINTED
C              IF K<0, CROSS SECTION IS PRINTED)
C       IN IS INPUT UNIT
C       IOUT IS OUTPUT UNIT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      ! dummy
      double precision, intent(inout) :: A(JJ,II)
      CHARACTER*24,     intent(inout) :: ANAME
      integer,          intent(in)    :: ii, jj, k
      integer,          intent(inout) :: in
      integer,          intent(in)    :: iout
      double precision, intent(inout), optional :: constantval
      logical,          intent(inout), optional :: constant
      integer,          intent(inout), optional :: iprnarg
      type(FileListType), pointer, optional :: Mf2005Files
      ! local
      integer :: istat
      double precision :: r, cnstnt
      character(len=1000) :: msg
      CHARACTER*20 FMTIN
      CHARACTER*200 CNTRL
      CHARACTER*16 TEXT
      CHARACTER*200 FNAME
      double precision :: PERTIM,TOTIM
      real :: pertims, totims
      real, allocatable, dimension(:,:) :: asingle
      type(FileType), pointer :: fil => null()
      DATA NUNOPN/99/
C     ------------------------------------------------------------------
C
      iprec = 0
C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
      READ(IN,'(A)') CNTRL
C
C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
      ICLOSE=0
      IFREE=1
      ICOL=1
      if (present(constant)) constant = .false.
      CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
      IF (CNTRL(ISTART:ISTOP).EQ.'CONSTANT') THEN
         LOCAT=0
         if (present(constant)) constant = .true.
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'INTERNAL') THEN
         LOCAT=IN
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'EXTERNAL') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,LOCAT,R,IOUT,IN)
      ELSE IF(CNTRL(ISTART:ISTOP).EQ.'OPEN/CLOSE') THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,0,N,R,IOUT,IN)
         FNAME=CNTRL(ISTART:ISTOP)
         LOCAT=NUNOPN
         WRITE(IOUT,15) LOCAT,FNAME
   15    FORMAT(1X,/1X,'OPENING FILE ON UNIT ',I4,':',/1X,A)
         ICLOSE=1
      ELSE
C
C2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
C2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
         IFREE=0
         READ(CNTRL,1,ERR=500) LOCAT,CNSTNT,FMTIN,IPRN
    1    FORMAT(I10,F10.0,A20,I10)
         if (present(constant) .and. locat==0) constant = .true.
         if (present(iprnarg)) then
           iprnarg = iprn
         endif
      END IF
C
C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
      IF(IFREE.NE.0) THEN
         CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,3,N,CNSTNT,IOUT,IN)
         IF(LOCAT.NE.0) THEN
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,1,N,R,IOUT,IN)
            FMTIN=CNTRL(ISTART:ISTOP)
            IF(ICLOSE.NE.0) THEN
               ! Data will be read from an OPEN/CLOSE file
               IF(FMTIN.EQ.'(BINARY)') THEN
                  OPEN(UNIT=LOCAT,FILE=FNAME,FORM=FORM,ACCESS=ACCESS,
     &                 ACTION=ACTION(1),iostat=istat)
                  if (istat /= 0) then
                    msg = 'Cannot open file: ' // trim(fname)
                    call store_error(msg)
                    call ustop()
                  endif
                  ! Determine precision (set iprec = 1 or 2) 
                  ! of a binary file by calling HEADPRECISION.
                  iprec = HEADPRECISION(locat)
               ELSE
                  OPEN(UNIT=LOCAT,FILE=FNAME,ACTION=ACTION(1),
     &                 iostat=istat)
                  if (istat /= 0) then
                    msg = 'Cannot open file: ' // trim(fname)
                    call store_error(msg)
                    call ustop()
                  endif
               END IF
            END IF
            IF(LOCAT.GT.0 .AND. FMTIN.EQ.'(BINARY)')then
                ! If precision of external (not OPEN/CLOSE) binary file has not been determined,
                ! determine precision (set iprec = 1 or 2) by calling HEADPRECISION.
                if (present(Mf2005Files)) then
                  fil => Mf2005Files%GetFileByUnit(locat)
                  if (associated(fil)) then
                    if (fil%IPrec > 0) then
                      iprec = fil%IPrec
                    else
                      iprec = HEADPRECISION(locat)
                      fil%IPrec = iprec
                    endif
                  endif
                endif
                LOCAT=-LOCAT
            endif
            CALL URWORD(CNTRL,ICOL,ISTART,ISTOP,2,IPRN,R,IOUT,IN)
            if (present(iprnarg)) then
               iprnarg = iprn
            endif
         END IF
      END IF
      !
      if (present(constantval)) then
        constantval = cnstnt
      endif
C
C4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
      IF(LOCAT.EQ.0) THEN
C
C4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
        if (present(constantval)) then
          constantval = cnstnt
        endif
        DO 80 I=1,II
        DO 80 J=1,JJ
   80   A(J,I)=CNSTNT
        IF(K.GT.0) WRITE(IOUT,2) ANAME,CNSTNT,K
    2   FORMAT(1X,/1X,A,' =',1P,G14.6,' FOR LAYER',I4)
        IF(K.LE.0) WRITE(IOUT,3) ANAME,CNSTNT
    3   FORMAT(1X,/1X,A,' =',1P,G14.6)
        RETURN
      ELSE IF(LOCAT.GT.0) THEN
C
C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
        IF(K.GT.0) THEN
           WRITE(IOUT,94) ANAME,K,LOCAT,FMTIN
   94      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
        ELSE IF(K.EQ.0) THEN
           WRITE(IOUT,95) ANAME,LOCAT,FMTIN
   95      FORMAT(1X,///11X,A,/
     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
        ELSE
           WRITE(IOUT,96) ANAME,LOCAT,FMTIN
   96      FORMAT(1X,///11X,A,' FOR CROSS SECTION',/
     1      1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
        END IF
        DO 100 I=1,II
        IF(FMTIN.EQ.'(FREE)') THEN
           READ(LOCAT,*) (A(J,I),J=1,JJ)
        ELSE
           READ(LOCAT,FMTIN) (A(J,I),J=1,JJ)
        END IF
  100   CONTINUE
      ELSE
C
C4C-----LOCAT<0; READ UNFORMATTED ARRAY VALUES.
        LOCAT=-LOCAT
        IF(K.GT.0) THEN
           WRITE(IOUT,201) ANAME,K,LOCAT
  201      FORMAT(1X,///11X,A,' FOR LAYER',I4,/
     1      1X,'READING BINARY ON UNIT ',I4)
        ELSE IF(K.EQ.0) THEN
           WRITE(IOUT,202) ANAME,LOCAT
  202      FORMAT(1X,///1X,A,/
     1      1X,'READING BINARY ON UNIT ',I4)
        ELSE
           WRITE(IOUT,203) ANAME,LOCAT
  203      FORMAT(1X,///1X,A,' FOR CROSS SECTION',/
     1      1X,'READING BINARY ON UNIT ',I4)
        END IF
        
        
        ! Ned todo: read binary file using single or double precision
        if (iprec == 1) then
           allocate(asingle(jj,ii))
           READ(LOCAT) KSTP,KPER,pertims,totims,TEXT,NCOL,NROW,ILAY
           PERTIM = real(pertims, kind(1.0d0))
           TOTIM = real(totims, kind(1.0d0))
           READ(LOCAT) asingle
           do i=1,ii
              do j=1,jj
                 a(j,i) = real(asingle(j,i), kind(1.0d0))
              enddo
           enddo
           deallocate(asingle)
        elseif (iprec == 2) then
           READ(LOCAT) KSTP,KPER,PERTIM,TOTIM,TEXT,NCOL,NROW,ILAY
           READ(LOCAT) A
        else
           ! Ned todo: should be an error.
        endif
      END IF
C
C5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
      IF(ICLOSE.NE.0) CLOSE(UNIT=LOCAT)
      ZERO=0.
      IF(CNSTNT.EQ.ZERO) GO TO 320
      DO 310 I=1,II
      DO 310 J=1,JJ
      A(J,I)=A(J,I)*CNSTNT
  310 CONTINUE
C
C6------IF PRINT CODE (IPRN) >0 OR =0 THEN PRINT ARRAY VALUES.
  320 IF(IPRN.GE.0) CALL ULAPRW(A,ANAME,0,0,JJ,II,0,IPRN,IOUT)
C
C7------RETURN
      RETURN
C
C8------CONTROL RECORD ERROR.
  500 IF(K.GT.0) THEN
         WRITE(IOUT,501) ANAME,K
  501    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,
     1     ' FOR LAYER',I4,':')
      ELSE
         WRITE(IOUT,502) ANAME
  502    FORMAT(1X,/1X,'ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
      END IF
      WRITE(IOUT,'(1X,A)') CNTRL
      CALL USTOPx(' ')
      END SUBROUTINE U2DREL
C*******************************************************************************
      SUBROUTINE UCOLNO(NLBL1,NLBL2,NSPACE,NCPL,NDIG,IOUT)
C     ******************************************************************
C     OUTPUT COLUMN NUMBERS ABOVE A MATRIX PRINTOUT
C        NLBL1 IS THE START COLUMN LABEL (NUMBER)
C        NLBL2 IS THE STOP COLUMN LABEL (NUMBER)
C        NSPACE IS NUMBER OF BLANK SPACES TO LEAVE AT START OF LINE
C        NCPL IS NUMBER OF COLUMN NUMBERS PER LINE
C        NDIG IS NUMBER OF CHARACTERS IN EACH COLUMN FIELD
C        IOUT IS OUTPUT CHANNEL
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*1 DOT,SPACE,DG,BF
      DIMENSION BF(130),DG(10)
C
      DATA DG(1),DG(2),DG(3),DG(4),DG(5),DG(6),DG(7),DG(8),DG(9),DG(10)/
     1         '0','1','2','3','4','5','6','7','8','9'/
      DATA DOT,SPACE/'.',' '/
C     ------------------------------------------------------------------
C
C1------CALCULATE # OF COLUMNS TO BE PRINTED (NLBL), WIDTH
C1------OF A LINE (NTOT), NUMBER OF LINES (NWRAP).
      WRITE(IOUT,1)
    1 FORMAT(1X)
      NLBL=NLBL2-NLBL1+1
      N=NLBL
      IF(NLBL.GT.NCPL) N=NCPL
      NTOT=NSPACE+N*NDIG
      IF(NTOT.GT.130) GO TO 50
      NWRAP=(NLBL-1)/NCPL + 1
      J1=NLBL1-NCPL
      J2=NLBL1-1
C
C2------BUILD AND PRINT EACH LINE
      DO 40 N=1,NWRAP
C
C3------CLEAR THE BUFFER (BF).
      DO 20 I=1,130
      BF(I)=SPACE
   20 CONTINUE
      NBF=NSPACE
C
C4------DETERMINE FIRST (J1) AND LAST (J2) COLUMN # FOR THIS LINE.
      J1=J1+NCPL
      J2=J2+NCPL
      IF(J2.GT.NLBL2) J2=NLBL2
C
C5------LOAD THE COLUMN #'S INTO THE BUFFER.
      DO 30 J=J1,J2
      NBF=NBF+NDIG
      I2=J/10
      I1=J-I2*10+1
      BF(NBF)=DG(I1)
      IF(I2.EQ.0) GO TO 30
      I3=I2/10
      I2=I2-I3*10+1
      BF(NBF-1)=DG(I2)
      IF(I3.EQ.0) GO TO 30
      I4=I3/10
      I3=I3-I4*10+1
      BF(NBF-2)=DG(I3)
      IF(I4.EQ.0) GO TO 30
      IF(I4.GT.9) THEN
C5A-----If more than 4 digits, use "X" for 4th digit.
         BF(NBF-3)='X'
      ELSE
         BF(NBF-3)=DG(I4+1)
      END IF
   30 CONTINUE
C
C6------PRINT THE CONTENTS OF THE BUFFER (I.E. PRINT THE LINE).
      WRITE(IOUT,31) (BF(I),I=1,NBF)
   31 FORMAT(1X,130A1)
C
   40 CONTINUE
C
C7------PRINT A LINE OF DOTS (FOR ESTHETIC PURPOSES ONLY).
   50 NTOT=NTOT
      IF(NTOT.GT.130) NTOT=130
      WRITE(IOUT,51) (DOT,I=1,NTOT)
   51 FORMAT(1X,130A1)
C
C8------RETURN
      RETURN
      END SUBROUTINE UCOLNO
C*******************************************************************************
      SUBROUTINE ULAPRS(BUF,TEXT,KSTP,KPER,NCOL,NROW,ILAY,IPRN,IOUT)
C     ******************************************************************
C     PRINT A 1 LAYER ARRAY IN STRIPS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DIMENSION BUF(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------MAKE SURE THE FORMAT CODE (IP OR IPRN) IS BETWEEN 1
C1------AND 21.
      IP=IPRN
      IF(IP.LT.1 .OR. IP.GT.21) IP=12
C
C2------DETERMINE THE NUMBER OF VALUES (NCAP) PRINTED ON ONE LINE.
      NCAP=10
      IF(IP.EQ.1) NCAP=11
      IF(IP.EQ.2) NCAP=9
      IF(IP.GT.2 .AND. IP.LT.7) NCAP=15
      IF(IP.GT.6 .AND. IP.LT.12) NCAP=20
      IF(IP.EQ.19) NCAP=5
      IF(IP.EQ.20) NCAP=6
      IF(IP.EQ.21) NCAP=7
C
C3------CALCULATE THE NUMBER OF STRIPS (NSTRIP).
      NCPF=129/NCAP
      IF(IP.GE.13 .AND. IP.LE.18) NCPF=7
      IF(IP.EQ.19) NCPF=13
      IF(IP.EQ.20) NCPF=12
      IF(IP.EQ.21) NCPF=10
      ISP=0
      IF(NCAP.GT.12 .OR. IP.GE.13) ISP=3
      NSTRIP=(NCOL-1)/NCAP + 1
      J1=1-NCAP
      J2=0
C
C4------LOOP THROUGH THE STRIPS.
      DO 2000 N=1,NSTRIP
C
C5------CALCULATE THE FIRST(J1) & THE LAST(J2) COLUMNS FOR THIS STRIP
      J1=J1+NCAP
      J2=J2+NCAP
      IF(J2.GT.NCOL) J2=NCOL
C
C6-------PRINT TITLE ON EACH STRIP DEPENDING ON ILAY
      IF(ILAY.GT.0) THEN
         WRITE(IOUT,1) TEXT,ILAY,KSTP,KPER
    1    FORMAT('1',/2X,A,' IN LAYER ',I3,' AT END OF TIME STEP ',I3,
     1     ' IN STRESS PERIOD ',I4/2X,75('-'))
      ELSE IF(ILAY.LT.0) THEN
         WRITE(IOUT,2) TEXT,KSTP,KPER
    2    FORMAT('1',/1X,A,' FOR CROSS SECTION AT END OF TIME STEP',I3,
     1     ' IN STRESS PERIOD ',I4/1X,79('-'))
      END IF
C
C7------PRINT COLUMN NUMBERS ABOVE THE STRIP
      CALL UCOLNO(J1,J2,ISP,NCAP,NCPF,IOUT)
C
C8------LOOP THROUGH THE ROWS PRINTING COLS J1 THRU J2 WITH FORMAT IP
      DO 1000 I=1,NROW
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     1      180,190,200,210), IP
C
C------------FORMAT 11G10.3
   10 WRITE(IOUT,11) I,(BUF(J,I),J=J1,J2)
   11 FORMAT(1X,I3,2X,1PG10.3,10(1X,G10.3))
      GO TO 1000
C
C------------FORMAT 9G13.6
   20 WRITE(IOUT,21) I,(BUF(J,I),J=J1,J2)
   21 FORMAT(1X,I3,2X,1PG13.6,8(1X,G13.6))
      GO TO 1000
C
C------------FORMAT 15F7.1
   30 WRITE(IOUT,31) I,(BUF(J,I),J=J1,J2)
   31 FORMAT(1X,I3,1X,15(1X,F7.1))
      GO TO 1000
C
C------------FORMAT 15F7.2
   40 WRITE(IOUT,41) I,(BUF(J,I),J=J1,J2)
   41 FORMAT(1X,I3,1X,15(1X,F7.2))
      GO TO 1000
C
C------------FORMAT 15F7.3
   50 WRITE(IOUT,51) I,(BUF(J,I),J=J1,J2)
   51 FORMAT(1X,I3,1X,15(1X,F7.3))
      GO TO 1000
C
C------------FORMAT 15F7.4
   60 WRITE(IOUT,61) I,(BUF(J,I),J=J1,J2)
   61 FORMAT(1X,I3,1X,15(1X,F7.4))
      GO TO 1000
C
C------------FORMAT 20F5.0
   70 WRITE(IOUT,71) I,(BUF(J,I),J=J1,J2)
   71 FORMAT(1X,I3,1X,20(1X,F5.0))
      GO TO 1000
C
C------------FORMAT 20F5.1
   80 WRITE(IOUT,81) I,(BUF(J,I),J=J1,J2)
   81 FORMAT(1X,I3,1X,20(1X,F5.1))
      GO TO 1000
C
C------------FORMAT 20F5.2
   90 WRITE(IOUT,91) I,(BUF(J,I),J=J1,J2)
   91 FORMAT(1X,I3,1X,20(1X,F5.2))
      GO TO 1000
C
C------------FORMAT 20F5.3
  100 WRITE(IOUT,101) I,(BUF(J,I),J=J1,J2)
  101 FORMAT(1X,I3,1X,20(1X,F5.3))
      GO TO 1000
C
C------------FORMAT 20F5.4
  110 WRITE(IOUT,111) I,(BUF(J,I),J=J1,J2)
  111 FORMAT(1X,I3,1X,20(1X,F5.4))
      GO TO 1000
C
C------------FORMAT 10G11.4
  120 WRITE(IOUT,121) I,(BUF(J,I),J=J1,J2)
  121 FORMAT(1X,I3,2X,1PG11.4,9(1X,G11.4))
      GO TO 1000
C
C------------FORMAT 10F6.0
  130 WRITE(IOUT,131) I,(BUF(J,I),J=J1,J2)
  131 FORMAT(1X,I3,1X,10(1X,F6.0))
      GO TO 1000
C
C------------FORMAT 10F6.1
  140 WRITE(IOUT,141) I,(BUF(J,I),J=J1,J2)
  141 FORMAT(1X,I3,1X,10(1X,F6.1))
      GO TO 1000
C
C------------FORMAT 10F6.2
  150 WRITE(IOUT,151) I,(BUF(J,I),J=J1,J2)
  151 FORMAT(1X,I3,1X,10(1X,F6.2))
      GO TO 1000
C
C------------FORMAT 10F6.3
  160 WRITE(IOUT,161) I,(BUF(J,I),J=J1,J2)
  161 FORMAT(1X,I3,1X,10(1X,F6.3))
      GO TO 1000
C
C------------FORMAT 10F6.4
  170 WRITE(IOUT,171) I,(BUF(J,I),J=J1,J2)
  171 FORMAT(1X,I3,1X,10(1X,F6.4))
      GO TO 1000
C
C------------FORMAT 10F6.5
  180 WRITE(IOUT,181) I,(BUF(J,I),J=J1,J2)
  181 FORMAT(1X,I3,1X,10(1X,F6.5))
C
C------------FORMAT 5G12.5
  190 WRITE(IOUT,191) I,(BUF(J,I),J=J1,J2)
  191 FORMAT(1X,I3,1X,1PG12.5,4(1X,G12.5))
C
C------------FORMAT 6G11.4
  200 WRITE(IOUT,201) I,(BUF(J,I),J=J1,J2)
  201 FORMAT(1X,I3,1X,1PG11.4,5(1X,G11.4))
C
C------------FORMAT 7G9.2
  210 WRITE(IOUT,211) I,(BUF(J,I),J=J1,J2)
  211 FORMAT(1X,I3,1X,1PG9.2,6(1X,G9.2))
C
 1000 CONTINUE
 2000 CONTINUE
C
C9------RETURN
      RETURN
      END SUBROUTINE ULAPRS
C*******************************************************************************
      SUBROUTINE ULAPRW(BUF,TEXT,KSTP,KPER,NCOL,NROW,ILAY,IPRN,IOUT)
C     ******************************************************************
C     PRINT 1 LAYER ARRAY
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      double precision :: BUF(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------PRINT A HEADER DEPENDING ON ILAY
      IF(ILAY.GT.0) THEN
         WRITE(IOUT,1) TEXT,ILAY,KSTP,KPER
    1    FORMAT('1',/2X,A,' IN LAYER ',I3,' AT END OF TIME STEP ',I3,
     1     ' IN STRESS PERIOD ',I4/2X,75('-'))
      ELSE IF(ILAY.LT.0) THEN
         WRITE(IOUT,2) TEXT,KSTP,KPER
    2    FORMAT('1',/1X,A,' FOR CROSS SECTION AT END OF TIME STEP',I3,
     1     ' IN STRESS PERIOD ',I4/1X,79('-'))
      END IF
C
C2------MAKE SURE THE FORMAT CODE (IP OR IPRN) IS
C2------BETWEEN 1 AND 21.
    5 IP=IPRN
      IF(IP.LT.1 .OR. IP.GT.21) IP=12
C
C3------CALL THE UTILITY MODULE UCOLNO TO PRINT COLUMN NUMBERS.
      IF(IP.EQ.1) CALL UCOLNO(1,NCOL,0,11,11,IOUT)
      IF(IP.EQ.2) CALL UCOLNO(1,NCOL,0,9,14,IOUT)
      IF(IP.GE.3 .AND. IP.LE.6) CALL UCOLNO(1,NCOL,3,15,8,IOUT)
      IF(IP.GE.7 .AND. IP.LE.11) CALL UCOLNO(1,NCOL,3,20,6,IOUT)
      IF(IP.EQ.12) CALL UCOLNO(1,NCOL,0,10,12,IOUT)
      IF(IP.GE.13 .AND. IP.LE.18) CALL UCOLNO(1,NCOL,3,10,7,IOUT)
      IF(IP.EQ.19) CALL UCOLNO(1,NCOL,0,5,13,IOUT)
      IF(IP.EQ.20) CALL UCOLNO(1,NCOL,0,6,12,IOUT)
      IF(IP.EQ.21) CALL UCOLNO(1,NCOL,0,7,10,IOUT)
C
C4------LOOP THROUGH THE ROWS PRINTING EACH ONE IN ITS ENTIRETY.
      DO 1000 I=1,NROW
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     1      180,190,200,210), IP
C
C------------ FORMAT 11G10.3
   10 WRITE(IOUT,11) I,(BUF(J,I),J=1,NCOL)
   11 FORMAT(1X,I3,2X,1PG10.3,10(1X,G10.3):/(5X,11(1X,G10.3)))
      GO TO 1000
C
C------------ FORMAT 9G13.6
   20 WRITE(IOUT,21) I,(BUF(J,I),J=1,NCOL)
   21 FORMAT(1X,I3,2X,1PG13.6,8(1X,G13.6):/(5X,9(1X,G13.6)))
      GO TO 1000
C
C------------ FORMAT 15F7.1
   30 WRITE(IOUT,31) I,(BUF(J,I),J=1,NCOL)
   31 FORMAT(1X,I3,1X,15(1X,F7.1):/(5X,15(1X,F7.1)))
      GO TO 1000
C
C------------ FORMAT 15F7.2
   40 WRITE(IOUT,41) I,(BUF(J,I),J=1,NCOL)
   41 FORMAT(1X,I3,1X,15(1X,F7.2):/(5X,15(1X,F7.2)))
      GO TO 1000
C
C------------ FORMAT 15F7.3
   50 WRITE(IOUT,51) I,(BUF(J,I),J=1,NCOL)
   51 FORMAT(1X,I3,1X,15(1X,F7.3):/(5X,15(1X,F7.3)))
      GO TO 1000
C
C------------ FORMAT 15F7.4
   60 WRITE(IOUT,61) I,(BUF(J,I),J=1,NCOL)
   61 FORMAT(1X,I3,1X,15(1X,F7.4):/(5X,15(1X,F7.4)))
      GO TO 1000
C
C------------ FORMAT 20F5.0
   70 WRITE(IOUT,71) I,(BUF(J,I),J=1,NCOL)
   71 FORMAT(1X,I3,1X,20(1X,F5.0):/(5X,20(1X,F5.0)))
      GO TO 1000
C
C------------ FORMAT 20F5.1
   80 WRITE(IOUT,81) I,(BUF(J,I),J=1,NCOL)
   81 FORMAT(1X,I3,1X,20(1X,F5.1):/(5X,20(1X,F5.1)))
      GO TO 1000
C
C------------ FORMAT 20F5.2
   90 WRITE(IOUT,91) I,(BUF(J,I),J=1,NCOL)
   91 FORMAT(1X,I3,1X,20(1X,F5.2):/(5X,20(1X,F5.2)))
      GO TO 1000
C
C------------ FORMAT 20F5.3
  100 WRITE(IOUT,101) I,(BUF(J,I),J=1,NCOL)
  101 FORMAT(1X,I3,1X,20(1X,F5.3):/(5X,20(1X,F5.3)))
      GO TO 1000
C
C------------ FORMAT 20F5.4
  110 WRITE(IOUT,111) I,(BUF(J,I),J=1,NCOL)
  111 FORMAT(1X,I3,1X,20(1X,F5.4):/(5X,20(1X,F5.4)))
      GO TO 1000
C
C------------ FORMAT 10G11.4
  120 WRITE(IOUT,121) I,(BUF(J,I),J=1,NCOL)
  121 FORMAT(1X,I3,2X,1PG11.4,9(1X,G11.4):/(5X,10(1X,G11.4)))
      GO TO 1000
C
C------------ FORMAT 10F6.0
  130 WRITE(IOUT,131) I,(BUF(J,I),J=1,NCOL)
  131 FORMAT(1X,I3,1X,10(1X,F6.0):/(5X,10(1X,F6.0)))
      GO TO 1000
C
C------------ FORMAT 10F6.1
  140 WRITE(IOUT,141) I,(BUF(J,I),J=1,NCOL)
  141 FORMAT(1X,I3,1X,10(1X,F6.1):/(5X,10(1X,F6.1)))
      GO TO 1000
C
C------------ FORMAT 10F6.2
  150 WRITE(IOUT,151) I,(BUF(J,I),J=1,NCOL)
  151 FORMAT(1X,I3,1X,10(1X,F6.2):/(5X,10(1X,F6.2)))
      GO TO 1000
C
C------------ FORMAT 10F6.3
  160 WRITE(IOUT,161) I,(BUF(J,I),J=1,NCOL)
  161 FORMAT(1X,I3,1X,10(1X,F6.3):/(5X,10(1X,F6.3)))
      GO TO 1000
C
C------------ FORMAT 10F6.4
  170 WRITE(IOUT,171) I,(BUF(J,I),J=1,NCOL)
  171 FORMAT(1X,I3,1X,10(1X,F6.4):/(5X,10(1X,F6.4)))
      GO TO 1000
C
C------------ FORMAT 10F6.5
  180 WRITE(IOUT,181) I,(BUF(J,I),J=1,NCOL)
  181 FORMAT(1X,I3,1X,10(1X,F6.5):/(5X,10(1X,F6.5)))
      GO TO 1000
C
C------------FORMAT 5G12.5
  190 WRITE(IOUT,191) I,(BUF(J,I),J=1,NCOL)
  191 FORMAT(1X,I3,2X,1PG12.5,4(1X,G12.5):/(5X,5(1X,G12.5)))
      GO TO 1000
C
C------------FORMAT 6G11.4
  200 WRITE(IOUT,201) I,(BUF(J,I),J=1,NCOL)
  201 FORMAT(1X,I3,2X,1PG11.4,5(1X,G11.4):/(5X,6(1X,G11.4)))
      GO TO 1000
C
C------------FORMAT 7G9.2
  210 WRITE(IOUT,211) I,(BUF(J,I),J=1,NCOL)
  211 FORMAT(1X,I3,2X,1PG9.2,6(1X,G9.2):/(5X,7(1X,G9.2)))
C
 1000 CONTINUE
C
C5------RETURN
      RETURN
      END SUBROUTINE ULAPRW
C*******************************************************************************
      SUBROUTINE ULAPRWC(A,NCOL,NROW,ILAY,IOUT,IPRN,ANAME)
C     ******************************************************************
C     WRITE A TWO-DIMENSIONAL REAL ARRAY.  IF THE ARRAY IS CONSTANT,
C     PRINT JUST THE CONSTANT VALUE.  IF THE ARRAY IS NOT CONSTANT, CALL
C     ULAPRW TO PRINT IT.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      double precision :: A(NCOL,NROW)
      CHARACTER*(*) ANAME
C     ------------------------------------------------------------------
C
C  Check to see if entire array is a constant.
      TMP=A(1,1)
      DO 300 I=1,NROW
      DO 300 J=1,NCOL
      IF(A(J,I).NE.TMP) GO TO 400
  300 CONTINUE
      IF(ILAY.GT.0) THEN
         WRITE(IOUT,302) ANAME,TMP,ILAY
  302    FORMAT(1X,/1X,A,' =',1P,G14.6,' FOR LAYER',I4)
      ELSE IF(ILAY.EQ.0) THEN
         WRITE(IOUT,303) ANAME,TMP
  303    FORMAT(1X,/1X,A,' =',1P,G14.6)
      ELSE
         WRITE(IOUT,304) ANAME,TMP
  304    FORMAT(1X,/1X,A,' =',1P,G14.6,' FOR CROSS SECTION')
      END IF
      RETURN
C
C  Print the array.
  400 IF(ILAY.GT.0) THEN
         WRITE(IOUT,494) ANAME,ILAY
  494    FORMAT(1X,//11X,A,' FOR LAYER',I4)
      ELSE IF(ILAY.EQ.0) THEN
         WRITE(IOUT,495) ANAME
  495    FORMAT(1X,//11X,A)
      ELSE
         WRITE(IOUT,496) ANAME
  496    FORMAT(1X,//11X,A,' FOR CROSS SECTION')
      END IF
      IF(IPRN.GE.0) CALL ULAPRW(A,ANAME,0,0,NCOL,NROW,0,IPRN,IOUT)
C
      RETURN
      END SUBROUTINE ULAPRWC
C*******************************************************************************
      SUBROUTINE ULASAV(BUF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                   NROW,ILAY,ICHN)
C     ******************************************************************
C     SAVE 1 LAYER ARRAY ON DISK
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      double precision :: BUF(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD CONTAINING IDENTIFYING
C1------INFORMATION.
      WRITE(ICHN) KSTP,KPER,PERTIM,TOTIM,TEXT,NCOL,NROW,ILAY
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING ARRAY VALUES
C2------THE ARRAY IS DIMENSIONED (NCOL,NROW)
      WRITE(ICHN) ((BUF(IC,IR),IC=1,NCOL),IR=1,NROW)
C
C3------RETURN
      RETURN
      END SUBROUTINE ULASAV
C*******************************************************************************
      SUBROUTINE ULASV2(BUFF,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                   NROW,ILAY,ICHN,FMTOUT,LBLSAV,IBOUND)
C     ******************************************************************
C     SAVE 1 LAYER ARRAY ON DISK USING FORMATTED OUTPUT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      double precision :: BUFF(NCOL,NROW)
      integer :: IBOUND(NCOL,NROW)
      CHARACTER*20 FMTOUT
C     ------------------------------------------------------------------
C
C1------WRITE A LABEL IF LBLSAV IS NOT 0.
      IF(LBLSAV.NE.0) WRITE(ICHN,5) KSTP,KPER,PERTIM,TOTIM,TEXT,NCOL,
     1                 NROW,ILAY,FMTOUT
5     FORMAT(1X,2I5,1P,2E15.6,1X,A,3I6,1X,A)
C
C2------WRITE THE ARRAY USING THE SPECIFIED FORMAT.
      DO 10 IR=1,NROW
      WRITE(ICHN,FMTOUT) (BUFF(IC,IR),IC=1,NCOL)
10    CONTINUE
C
C3------RETURN
      RETURN
      END SUBROUTINE ULASV2
C*******************************************************************************
      SUBROUTINE ULASV3(IDATA,TEXT,KSTP,KPER,PERTIM,TOTIM,NCOL,
     1                   NROW,ILAY,ICHN,FMTOUT,LBLSAV)
C     ******************************************************************
C     SAVE 2-D (LAYER) INTEGER ARRAY ON DISK USING FORMATTED OUTPUT
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DIMENSION IDATA(NCOL,NROW)
      CHARACTER*20 FMTOUT
C     ------------------------------------------------------------------
C
C1------WRITE A LABEL IF LBLSAV IS NOT 0.
      IF(LBLSAV.NE.0) WRITE(ICHN,5) KSTP,KPER,PERTIM,TOTIM,TEXT,NCOL,
     1                 NROW,ILAY,FMTOUT
5     FORMAT(1X,2I5,1P,2E15.6,1X,A,3I6,1X,A)
C
C2------WRITE THE ARRAY USING THE SPECIFIED FORMAT.
      DO 10 IR=1,NROW
      WRITE(ICHN,FMTOUT) (IDATA(IC,IR),IC=1,NCOL)
10    CONTINUE
C
C3------RETURN.
      RETURN
      END SUBROUTINE ULASV3
C*******************************************************************************
      SUBROUTINE UBUDSV(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,IOUT)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      double precision :: BUFF(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE AN UNFORMATTED RECORD IDENTIFYING DATA.
      WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBUDSV SAVING "',A16,'" ON UNIT',I3,
     1     ' AT TIME STEP',I3,', STRESS PERIOD ',I4)
      WRITE(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,NLAY
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      WRITE(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END SUBROUTINE UBUDSV
C*******************************************************************************
      SUBROUTINE UBDSV1(KSTP,KPER,TEXT,IBDCHN,BUFF,NCOL,NROW,NLAY,IOUT,
     1          DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW AS A 3-D
C     ARRAY WITH EXTRA RECORD TO INDICATE DELT, PERTIM, AND TOTIM.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      double precision :: BUFF(NCOL,NROW,NLAY)
      integer :: IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE TWO UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV1 SAVING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I3,', STRESS PERIOD',I4)
      WRITE(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,-NLAY
      WRITE(IBDCHN) 1,DELT,PERTIM,TOTIM
C
C2------WRITE AN UNFORMATTED RECORD CONTAINING VALUES FOR
C2------EACH CELL IN THE GRID.
      WRITE(IBDCHN) BUFF
C
C3------RETURN
      RETURN
      END SUBROUTINE UBDSV1
C*******************************************************************************
      SUBROUTINE UBDSV2(KSTP,KPER,TEXT,IBDCHN,NCOL,NROW,NLAY,
     1          NLIST,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     WRITE HEADER RECORDS FOR CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT
C     OF FLOW USING A LIST STRUCTURE.  EACH ITEM IN THE LIST IS WRITTEN
C     BY MODULE UBDSVA
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      DIMENSION IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE THREE UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV2 SAVING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I3,', STRESS PERIOD',I4)
      WRITE(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,-NLAY
      WRITE(IBDCHN) 2,DELT,PERTIM,TOTIM
      WRITE(IBDCHN) NLIST
C
C2------RETURN
      RETURN
      END SUBROUTINE UBDSV2
C*******************************************************************************
      SUBROUTINE UBDSVA(IBDCHN,NCOL,NROW,J,I,K,Q,IBOUND,NLAY)
C     ******************************************************************
C     WRITE ONE VALUE OF CELL-BY-CELL FLOW USING A LIST STRUCTURE.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------CALCULATE CELL NUMBER
      ICRL= (K-1)*NROW*NCOL + (I-1)*NCOL + J
C
C2------WRITE CELL NUMBER AND FLOW RATE
      WRITE(IBDCHN) ICRL,Q
C
C3------RETURN
      RETURN
      END SUBROUTINE UBDSVA
C*******************************************************************************
      SUBROUTINE UBDSV3(KSTP,KPER,TEXT,IBDCHN,BUFF,IBUFF,NOPT,
     1              NCOL,NROW,NLAY,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     RECORD CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT OF FLOW AS A 2-D
C     ARRAY OF FLOW VALUES AND OPTIONALLY A 2-D ARRAY OF LAYER NUMBERS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT
      double precision :: BUFF(NCOL,NROW,NLAY)
      integer :: IBUFF(NCOL,NROW),
     1          IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE TWO UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV3 SAVING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I3,', STRESS PERIOD',I4)
      WRITE(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,-NLAY
      IMETH=3
      IF(NOPT.EQ.1) IMETH=4
      WRITE(IBDCHN) IMETH,DELT,PERTIM,TOTIM
C
C2------WRITE DATA AS ONE OR TWO UNFORMATTED RECORDS CONTAINING ONE
C2------VALUE PER LAYER.
      IF(NOPT.EQ.1) THEN
C2A-----WRITE ONE RECORD WHEN NOPT IS 1.  THE VALUES ARE FLOW VALUES
C2A-----FOR LAYER 1.
         WRITE(IBDCHN) ((BUFF(J,I,1),J=1,NCOL),I=1,NROW)
      ELSE
C2B-----WRITE TWO RECORDS WHEN NOPT IS NOT 1.  FIRST RECORD CONTAINS
C2B-----LAYER NUMBERS;  SECOND RECORD CONTAINS FLOW VALUES.
         WRITE(IBDCHN) ((IBUFF(J,I),J=1,NCOL),I=1,NROW)
         WRITE(IBDCHN) ((BUFF(J,I,IBUFF(J,I)),J=1,NCOL),I=1,NROW)
      END IF
C
C3------RETURN
      RETURN
      END SUBROUTINE UBDSV3
C*******************************************************************************
      SUBROUTINE UBDSV4(KSTP,KPER,TEXT,NAUX,AUXTXT,IBDCHN,
     1          NCOL,NROW,NLAY,NLIST,IOUT,DELT,PERTIM,TOTIM,IBOUND)
C     ******************************************************************
C     WRITE HEADER RECORDS FOR CELL-BY-CELL FLOW TERMS FOR ONE COMPONENT
C     OF FLOW PLUS AUXILIARY DATA USING A LIST STRUCTURE.  EACH ITEM IN
C     THE LIST IS WRITTEN BY MODULE UBDSVB
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*16 TEXT,AUXTXT(*)
      DIMENSION IBOUND(NCOL,NROW,NLAY)
C     ------------------------------------------------------------------
C
C1------WRITE UNFORMATTED RECORDS IDENTIFYING DATA.
      IF(IOUT.GT.0) WRITE(IOUT,1) TEXT,IBDCHN,KSTP,KPER
    1 FORMAT(1X,'UBDSV4 SAVING "',A16,'" ON UNIT',I4,
     1     ' AT TIME STEP',I3,', STRESS PERIOD',I4)
      WRITE(IBDCHN) KSTP,KPER,TEXT,NCOL,NROW,-NLAY
      WRITE(IBDCHN) 5,DELT,PERTIM,TOTIM
      WRITE(IBDCHN) NAUX+1
      IF(NAUX.GT.0) WRITE(IBDCHN) (AUXTXT(N),N=1,NAUX)
      WRITE(IBDCHN) NLIST
C
C2------RETURN
      RETURN
      END SUBROUTINE UBDSV4
C*******************************************************************************
      SUBROUTINE UBDSVB(IBDCHN,NCOL,NROW,J,I,K,Q,VAL,NVL,NAUX,LAUX,
     1                  IBOUND,NLAY)
C     ******************************************************************
C     WRITE ONE VALUE OF CELL-BY-CELL FLOW PLUS AUXILIARY DATA USING
C     A LIST STRUCTURE.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION IBOUND(NCOL,NROW,NLAY),VAL(NVL)
C     ------------------------------------------------------------------
C
C1------CALCULATE CELL NUMBER
      ICRL= (K-1)*NROW*NCOL + (I-1)*NCOL + J
C
C2------WRITE CELL NUMBER AND FLOW RATE
      IF(NAUX.GT.0) THEN
         N2=LAUX+NAUX-1
         WRITE(IBDCHN) ICRL,Q,(VAL(N),N=LAUX,N2)
      ELSE
         WRITE(IBDCHN) ICRL,Q
      END IF
C
C3------RETURN
      RETURN
      END SUBROUTINE UBDSVB
C*******************************************************************************
      SUBROUTINE UMESPR(TEXT1,TEXT2,IOUT)
C     ******************************************************************
C     PRINT A LINE CONSISTING OF TWO TEXT VARIABLES.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*(*) TEXT1,TEXT2
C     ------------------------------------------------------------------
      WRITE(IOUT,*)
      WRITE(IOUT,'(1X,2A)') TEXT1,TEXT2
C
      RETURN
      END SUBROUTINE UMESPR
C*******************************************************************************
      SUBROUTINE USTOPx(STOPMESS)
C     ******************************************************************
C     STOP PROGRAM, WITH OPTION TO PRINT MESSAGE BEFORE STOPPING
C     ******************************************************************
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER STOPMESS*(*)
C     ------------------------------------------------------------------
   10 FORMAT(1X,A)

      IF (STOPMESS.NE.' ') THEN
        WRITE(*,10) STOPMESS
      ENDIF
      STOP

      END SUBROUTINE USTOPx
  
      function build_obsname(onamebase,i,nchar) result(oname)
        ! Return a string of form: NAMEiiiiiii
        ! Where: NAME is trim(onamebase), and
        !        iiiiiii is integer i including leading zeros 
        !        nchar is number of characters to fill in writing iiiiiii
        implicit none
        ! dummy
        character(len=*), intent(in) :: onamebase
        integer,          intent(in) :: i
        integer,          intent(in) :: nchar
        character(len=LENOBSNAME)    :: oname
        ! local
        character(len=LENOBSNAME-1) :: cnum
        character(len=8)            :: numfmt
        ! format
   10   format('(i',i0,'.',i0,')')
        !
        ! Define format for writing i
        write(numfmt,10)nchar,nchar
        !
        ! Make observation name from NAME part and integer part.
        write(cnum,numfmt)i
        oname = trim(onamebase) // cnum
        !
        return
      end function build_obsname

      subroutine assign_ncharsizes_flow(ncharid)
        implicit none
        ! dummy
        integer, intent(inout) :: ncharid
        ! local
        integer :: lenmax=20, numres=5
        !
        ! Assign sizes for observations other than head
        ncharid = lenmax - numres
        !
        return
      end subroutine assign_ncharsizes_flow

      end module utl7module
