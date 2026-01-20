module ArrayReadersMF5Module
  
  use ConstantsModule,   only: LINELENGTH, LENBIGLINE, LENBOUNDNAME, &
                               NAMEDBOUNDFLAG, LINELENGTH, DZERO
  use InputOutputModule, only: openfile, u9rdcom, urword, ucolno, ulaprw
  use KindModule,        only: DP, I4B
  use OpenSpecModule,    only: ACCESS, FORM
  use SimPHMFModule,     only: store_error, ustop, store_error_unit

  implicit none

  private
  public :: ReadArray
  public :: u1ddbl, u1dint, u2ddbl, u2dint, u3ddbl, u3dint
  
  interface ReadArray
    module procedure u1ddbl, u1dint, u2ddbl, u2dint, u3ddbl, u3dint, &
                     u1ddbl_3darg, u1dint_3darg, u2dint_1darg, &
                     u3dint_1darg, u2ddbl_1darg, u3ddbl_1darg
  end interface ReadArray

! -- Integer array readers
! subroutine u1dint(ia,aname,jj,in,iout)        ! ia is a 1-D array
! subroutine u1dint_3darg(ia,aname,jj,in,iout)  ! ia is a 3-D array to be read as a 1-D data set
! subroutine u2dint(ia,aname,ii,jj,k,in,iout)              ! ia is a 2-D array to be read as a 1-D data set
! subroutine u2dint_1darg(ia, aname, ii, jj, k, in, iout)  ! ia is a 1-D array to be read as a 1-D data set
! subroutine u3dint(iin,iout,nlay,nrow,ncol,neq,ival,cval)

! -- Real array readers
! subroutine u1ddbl(a,aname,jj,in,iout)         ! a is a 1-D array
! subroutine u1ddbl_3darg(a,aname,jj,in,iout)   ! a is a 3-D array to be read as a 1-D data set
! subroutine u2ddbl(a,aname,ii,jj,k,in,iout)        ! a is a 2-D array to be read as a 1-D data set
! subroutine u2ddbl_1darg(a,aname,ii,jj,k,in,iout)  ! a is a 1-D array to be read as a 1-D data set
! subroutine u3ddbl(iin,iout,nlay,nrow,ncol,neq,rval,cval)       ! rval is a 3-D array to be read as a 1-D data set
! subroutine u3ddbl_1darg(iin,iout,nlay,nrow,ncol,neq,rval,cval) ! rval is a 1-D array to be read as a 1-D data set

contains

  subroutine u1ddbl_3darg(a,aname,jj,in,iout)
    ! Wrapper for u1ddbl to allow array argument 
    ! of ReadArray interface to be 3-D.
    ! dummy
    integer,                     intent(in)    :: jj
    real(DP), dimension(jj,1,1), intent(inout) :: a
    character(len=24),           intent(inout) :: aname
    integer,                     intent(in)    :: in, iout
    !
    call u1ddbl(a,aname,jj,in,iout)
    !
    return
  end subroutine u1ddbl_3darg

  subroutine u1ddbl(a,aname,jj,in,iout)
!   ******************************************************************
!   ROUTINE TO INPUT 1-D REAL DATA MATRICES
!     A IS ARRAY TO INPUT
!     ANAME IS 24 CHARACTER DESCRIPTION OF A
!     JJ IS NO. OF ELEMENTS
!     IN IS INPUT UNIT
!     IOUT IS OUTPUT UNIT
!   ******************************************************************
!
!   SPECIFICATIONS:
!   ------------------------------------------------------------------
    ! dummy
    integer, intent(in)                    :: jj
    real(DP), dimension(jj), intent(inout) :: a
    character(len=24), intent(in)          :: aname
    integer, intent(in)                    :: in, iout
    ! local
    integer :: iclose, icol, ierr, ifree, iprn, istart, istop, j, locat, n
    real(DP) :: r, cnstnt
    character(len=20) :: fmtin
    character (len=:), allocatable :: cntrl
    character(len=200) :: fname
    character(len=linelength) :: ermsg
    integer, parameter :: nunopn=99
    ! formats
    1 format(i10,f10.0,a20,i10)
    3 format(1x,/1x,a,' =',1p,g14.6)
    5 format(1x,///11x,a,/ &
             1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A20)
    502 format(' ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
    1001 format((1x,1pg12.5,9(1x,g12.5)))
    1002 format((1x,1pg12.5,4(1x,g12.5)))
!   ------------------------------------------------------------------
!
!C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
    call u9rdcom(in,iout,cntrl,ierr)
!
!C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
    iclose = 0
    ifree = 1
    icol = 1
    call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
    if (cntrl(istart:istop).eq.'CONSTANT') then
      locat = 0
    elseif (cntrl(istart:istop).eq.'INTERNAL') then
      locat = in
    elseif (cntrl(istart:istop).eq.'EXTERNAL') then
      write(ermsg, *) 'Error. EXTERNAL arrays not supported in MODFLOW 6.'
      call store_error(ermsg)
      write(ermsg, *) 'Use OPEN/CLOSE instead.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
    elseif (cntrl(istart:istop).eq.'OPEN/CLOSE') then
      call urword(cntrl,icol,istart,istop,0,n,r,iout,in)
      fname = cntrl(istart:istop)
      locat = nunopn
      call openfile(locat, iout, fname, 'OPEN/CLOSE')
      iclose=1
    else
      write(ermsg, *) 'ERROR READING CONTROL RECORD FOR ' //                 &
                       trim(adjustl(aname))
      call store_error(ermsg)
      call store_error(trim(adjustl(cntrl)))
      write(ermsg, *) 'Use CONSTANT, INTERNAL, or OPEN/CLOSE.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
!
!C2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!C2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
      ifree = 0
      read(cntrl,1,err=500) locat,cnstnt,fmtin,iprn
    endif
!
!C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
    if (ifree.ne.0) then
      call urword(cntrl,icol,istart,istop,3,n,cnstnt,iout,in)
      if (locat.gt.0) then
        call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
        fmtin = cntrl(istart:istop)
        call urword(cntrl,icol,istart,istop,2,iprn,r,iout,in)
      endif
    endif
!
!C4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
    if (locat.le.0) then
!
!C4A-----LOCAT <0 OR =0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
      do j=1,jj
        a(j) = cnstnt
      enddo
      if (iout>0) write(iout,3) aname,cnstnt
      return
    endif
!
!C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
    if (iout>0) write(iout,5) aname,locat,fmtin
    if (fmtin.eq.'(FREE)') then
      read(locat,*) (a(j),j=1,jj)
    else
      read(locat,fmtin) (a(j),j=1,jj)
    endif
    if (iclose.ne.0) close(unit=locat)
!
!C5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
    if (cnstnt.eq.DZERO) go to 120
    do j=1,jj
      a(j)=a(j)*cnstnt
    enddo
!
!C6------IF PRINT CODE (IPRN) =0 OR >0 THEN PRINT ARRAY VALUES.
    120 continue
    if (iout>0) then
      if (iprn.eq.0) then
        write(iout,1001) (a(j),j=1,jj)
      elseif (iprn.gt.0) then
        write(iout,1002) (a(j),j=1,jj)
      endif
    endif
!
!C7------RETURN
    return
!C
!C8------CONTROL RECORD ERROR.
    500 continue
    write(ermsg,502) aname
    call store_error(ermsg)
    write(ermsg,'(1x,a)') cntrl
    call store_error(ermsg)
    call store_error_unit(in)
    call ustop()
  end subroutine u1ddbl

  subroutine u1dint_3darg(a,aname,jj,in,iout)
    ! Wrapper for u1dint to allow array argument 
    ! of ReadArray interface to be 3-D.
    ! dummy
    integer,                         intent(in)    :: jj
    integer(I4B), dimension(jj,1,1), intent(inout) :: a
    character(len=24),               intent(inout) :: aname
    integer,                         intent(in)    :: in, iout
    !
    call u1dint(a,aname,jj,in,iout)
    !
    return
  end subroutine u1dint_3darg

  subroutine u1dint(ia,aname,jj,in,iout)
!   ******************************************************************
!   ROUTINE TO INPUT 1-D REAL DATA MATRICES
!     A IS ARRAY TO INPUT
!     ANAME IS 24 CHARACTER DESCRIPTION OF A
!     JJ IS NO. OF ELEMENTS
!     IN IS INPUT UNIT
!     IOUT IS OUTPUT UNIT
!   ******************************************************************
!
!      SPECIFICATIONS:
!   ------------------------------------------------------------------
    ! dummy
    integer, intent(in) :: jj
    integer, dimension(jj), intent(inout) :: ia
    character(len=24) :: aname
    integer, intent(in) :: in, iout
    ! local
    integer :: iclose, icnstnt, icol, ierr, ifree, iprn, istart, istop, j, n, &
               locat
    real(DP) :: r
    character(len=20) :: fmtin
    character (len=:), allocatable :: cntrl
    character(len=200) :: fname
    character(len=linelength) :: ermsg
    integer, parameter :: nunopn=99
    ! formats
    1 format(i10,i10.0,a20,i10)
    3 format(1x,/1x,a,' =',1p,i10)
    5 format(1X,///11X,A,/ &
            1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A20)
    502 format(' ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
    1001 format(20(1x,i9))
    1002 format(8(1x,i9))
!   ------------------------------------------------------------------
!
!C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
    call u9rdcom(in,iout,cntrl,ierr)
!C
!C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
    iclose=0
    ifree=1
    icol=1
    call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
    if (cntrl(istart:istop).eq.'CONSTANT') then
      locat=0
    elseif (cntrl(istart:istop).eq.'INTERNAL') then
      locat=in
    elseif (cntrl(istart:istop).eq.'EXTERNAL') then
      write(ermsg, *) 'Error. EXTERNAL arrays not supported in MODFLOW 6.'
      call store_error(ermsg)
      write(ermsg, *) 'Use OPEN/CLOSE instead.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
    elseif (cntrl(istart:istop).eq.'OPEN/CLOSE') then
      call urword(cntrl,icol,istart,istop,0,n,r,iout,in)
      fname=cntrl(istart:istop)
      locat=nunopn
      call openfile(locat, iout, fname, 'OPEN/CLOSE')
      iclose=1
    else
      write(ermsg, *) 'ERROR READING CONTROL RECORD FOR ' //                 &
                       trim(adjustl(aname))
      call store_error(ermsg)
      call store_error(trim(adjustl(cntrl)))
      write(ermsg, *) 'Use CONSTANT, INTERNAL, or OPEN/CLOSE.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
!C
!C2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!C2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
      ifree=0
      read(cntrl,1,err=500) locat,icnstnt,fmtin,iprn
    endif
!C
!C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
    if (ifree.ne.0) then
      call urword(cntrl,icol,istart,istop,2,icnstnt,r,iout,in)
      if (locat.gt.0) then
        call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
        fmtin=cntrl(istart:istop)
        call urword(cntrl,icol,istart,istop,2,iprn,r,iout,in)
      endif
    endif
!c
!c4------test locat to see how to define array values.
    if (locat.le.0) then
!C
!C4A-----LOCAT <0 OR =0; SET ALL ARRAY VALUES EQUAL TO ICNSTNT. RETURN.
      do j=1,jj
        ia(j)=icnstnt
      enddo
      if (iout>0) write(iout,3) aname,icnstnt
      return
    endif
!C
!C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
    if (iout>0) write(iout,5) aname,locat,fmtin
    if (fmtin.eq.'(FREE)') then
      read(locat,*) (ia(j),j=1,jj)
    else
      read(locat,fmtin) (ia(j),j=1,jj)
    endif
    if (iclose.ne.0) close(unit=locat)
!C
!C5------IF ICNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY ICNSTNT.
    if (icnstnt.ne.0) then
      do j=1,jj
        ia(j)=ia(j)*icnstnt
      enddo
    endif
!C
!C6------IF PRINT CODE (IPRN) =0 OR >0 THEN PRINT ARRAY VALUES.
    if (iout>0) then
      if (iprn.eq.0) then
        write(iout,1001) (ia(j),j=1,jj)
      elseif (iprn.gt.0) then
        write(iout,1002) (ia(j),j=1,jj)
      endif
    endif
!C
!C7------RETURN
    return
!C
!C8------CONTROL RECORD ERROR.
    500 continue
    write(ermsg,502) aname
    call store_error(ermsg)
    write(ermsg,'(1x,a)') cntrl
    call store_error(ermsg)
    call store_error_unit(in)
    call ustop()
  end subroutine u1dint

  subroutine u2ddbl_1darg(a,aname,ii,jj,k,in,iout)
    ! dummy
    real(DP),dimension(jj*ii), intent(inout) :: a
    character(len=24), intent(in) :: aname
    integer, intent(in) :: ii, jj, k, in, iout
    !
    call u2ddbl(a, aname, ii, jj, k, in, iout)
    !
    return
  end subroutine u2ddbl_1darg

  subroutine u2ddbl(a,aname,ii,jj,k,in,iout)
! ******************************************************************
! ROUTINE TO INPUT 2-D DOUBLE PRECISION DATA MATRICES
!   A IS ARRAY TO INPUT
!   ANAME IS 24 CHARACTER DESCRIPTION OF A
!   II IS NO. OF ROWS
!   JJ IS NO. OF COLS
!   K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --)
!          IF K=0, NO LAYER IS PRINTED
!          IF K<0, CROSS SECTION IS PRINTED)
!   IN IS INPUT UNIT
!   IOUT IS OUTPUT UNIT
! ******************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------
    ! dummy
    real(DP),dimension(jj,ii), intent(inout) :: a
    character(len=24), intent(in) :: aname
    integer, intent(in) :: ii, jj, k, in, iout
    ! local
    integer :: i, j, iclose, icol, ierr, ifree, ilay, iprn, istart, istop, &
               kper, kstp, locat, n, ncol, nrow
    real(DP) :: r, cnstnt, pertim, totim
    character(len=20) :: fmtin
    character (len=:), allocatable :: cntrl
    character(len=16) :: text
    character(len=200) :: fname
    character(len=20) :: ftype
    character(len=linelength) :: ermsg
    integer, parameter :: nunopn=99
    ! formats
    2 format(1x,/1x,a,' =',1p,g14.6,' FOR LAYER',i4)
    3 format(1x,/1x,a,' =',1p,g14.6)
    94 format(1x,///11x,a,' FOR LAYER',I4,/ &
              1X,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
    95 format(1x,///11x,a,/ &
              1x,'READING ON UNIT ',I4,' WITH FORMAT: ',A)
    96 format(1x,///11x,a,' FOR CROSS SECTION',/ &
              1X,'READING ON UNIT ',I4,' WITH FORMAT: ',a)
    201 format(1x,///11x,a,' FOR LAYER',i4,/ &
               1X,'READING BINARY ON UNIT ',i4)
    202 format(1x,///1x,a,/ &
               1x,'READING BINARY ON UNIT ',I4)
    203 format(1x,///1x,a,' FOR CROSS SECTION',/ &
               1X,'READING BINARY ON UNIT ',i4)
    501 format(' ERROR READING ARRAY CONTROL RECORD FOR ',A, &
               ' FOR LAYER',i4,':')
    502 format(' ERROR READING ARRAY CONTROL RECORD FOR ',A,':')
!     ------------------------------------------------------------------
!
!C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
    call u9rdcom(in,iout,cntrl,ierr)
!
!C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
    iclose=0
    ifree=1
    icol=1
    call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
    if (cntrl(istart:istop).eq.'CONSTANT') then
      locat=0
    elseif (cntrl(istart:istop).eq.'INTERNAL') then
      locat=in
    elseif (cntrl(istart:istop).eq.'EXTERNAL') then
      write(ermsg, *) 'Error. EXTERNAL arrays not supported in MODFLOW 6.'
      call store_error(ermsg)
      write(ermsg, *) 'Use OPEN/CLOSE instead.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
    elseif (cntrl(istart:istop).eq.'OPEN/CLOSE') then
      call urword(cntrl,icol,istart,istop,0,n,r,iout,in)
      fname=cntrl(istart:istop)
      locat=nunopn
      iclose=1
      ftype = 'OPEN/CLOSE'
    else
      write(ermsg, *) 'ERROR READING CONTROL RECORD FOR ' //                 &
                       trim(adjustl(aname))
      call store_error(ermsg)
      call store_error(trim(adjustl(cntrl)))
      write(ermsg, *) 'Use CONSTANT, INTERNAL, or OPEN/CLOSE.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
    endif
!
!C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
    if (ifree.ne.0) then
      call urword(cntrl,icol,istart,istop,3,n,cnstnt,iout,in)
      if (locat.ne.0) then
        call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
        fmtin=cntrl(istart:istop)
        if (iclose.ne.0) then
          if (fmtin.eq.'(BINARY)') then
            call openfile(locat, iout, fname, ftype, fmtarg_opt=FORM,    &
                          accarg_opt=ACCESS)
          else
            call openfile(locat, iout, fname, ftype)
          endif
        endif
        if (locat.gt.0 .and. fmtin.eq.'(BINARY)') locat=-locat
        call urword(cntrl,icol,istart,istop,2,iprn,r,iout,in)
      endif
    endif
!
!C4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
    if (locat.eq.0) then
!
!C4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO CNSTNT. RETURN.
      do i=1,ii
        do j=1,jj
          a(j,i)=cnstnt
        enddo
      enddo
      if (k.gt.0 .and. iout>0) write(iout,2) aname,cnstnt,k
      if (k.le.0 .and. iout>0) write(iout,3) aname,cnstnt
      return
    elseif (locat.gt.0) then
!
!C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
        if (iout>0) then
          if (k.gt.0) then
            write(iout,94) aname,k,locat,fmtin
          elseif (k.eq.0) then
            write(iout,95) aname,locat,fmtin
          else
            write(iout,96) aname,locat,fmtin
          endif
        endif
        do i=1,ii
          if (fmtin.eq.'(FREE)') then
            read(locat,*) (a(j,i),j=1,jj)
          else
            read(locat,fmtin) (a(j,i),j=1,jj)
          endif
        enddo
      else
!
!C4C-----LOCAT<0; READ UNFORMATTED ARRAY VALUES.
        locat=-locat
        if (iout>0) then
          if (k.gt.0) then
            write(iout,201) aname,k,locat
          elseif (k.eq.0) then
            write(iout,202) aname,locat
          else
            write(iout,203) aname,locat
          endif
        endif
        read(locat) kstp,kper,pertim,totim,text,ncol,nrow,ilay
        read(locat) a
      endif
!
!C5------IF CNSTNT NOT ZERO THEN MULTIPLY ARRAY VALUES BY CNSTNT.
      if (iclose.ne.0) close(unit=locat)
      if (cnstnt.ne.DZERO) then
        do i=1,ii
          do j=1,jj
            a(j,i)=a(j,i)*cnstnt
          enddo
        enddo
      endif
!
!C6------IF PRINT CODE (IPRN) >0 OR =0 THEN PRINT ARRAY VALUES.
!  320 continue
      if (iprn.ge.0) call ulaprw(a,aname,0,0,jj,ii,0,iprn,iout)
!
!C7------RETURN
      return
!
!C8------CONTROL RECORD ERROR.
      if (k.gt.0) then
        write(ermsg,501) aname,k
      else
        write(ermsg,502) aname
      endif
      call store_error(ermsg)
      write(ermsg,'(1x,a)') cntrl
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
  end subroutine u2ddbl
  
  subroutine u2dint_1darg(ia, aname, ii, jj, k, in, iout)
    ! dummy
    character(len=24),         intent(in) :: aname
    integer, dimension(jj*ii), intent(inout) :: ia
    integer,                   intent(in) :: ii, jj, k, in, iout
    !
    call u2dint(ia, aname, ii, jj, k, in, iout)
    !
    return
  end subroutine u2dint_1darg

  subroutine u2dint(ia,aname,ii,jj,k,in,iout)
! ******************************************************************
! ROUTINE TO INPUT 2-D INTEGER DATA MATRICES
!   IA IS ARRAY TO INPUT
!   ANAME IS 24 CHARACTER DESCRIPTION OF IA
!   II IS NO. OF ROWS
!   JJ IS NO. OF COLS
!   K IS LAYER NO. (USED WITH NAME TO TITLE PRINTOUT --
!          IF K=0, NO LAYER IS PRINTED
!          IF K<0, CROSS SECTION IS PRINTED)
!   IN IS INPUT UNIT
!   IOUT IS OUTPUT UNIT
! ******************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------
    ! dummy
    character(len=24), intent(in) :: aname
    integer, dimension(jj,ii), intent(inout) :: ia
    integer, intent(in) :: ii, jj, k, in, iout
    ! local
    integer :: i, j, iclose, icol, iconst, ierr, ifree, iprn, istart, &
               istop, locat, n
    real(DP) :: r
    character(len=20) :: fmtin
    character (len=:), allocatable :: cntrl
    character(len=200) :: fname
    character(len=20) :: ftype
    character(len=LINELENGTH) :: ermsg
    integer, parameter :: NUNOPN=99
    ! formats
    1 format(i10,i10,a20,i10)
    82 format(1x,/1x,a,' =',i15,' FOR LAYER',i4)
    83 format(1x,/1x,a,' =',i15)
    94 format(1x,///11x,a,' FOR LAYER',i4,/ &
              1x,'READING ON UNIT ',i4,' WITH FORMAT: ',a)
    95 format(1x,///11x,a,/ &
              1x,'READING ON UNIT ',I4,' WITH FORMAT: ',a)
    96 format(1x,///11x,a,' FOR CROSS SECTION',/ &
              1X,'READING ON UNIT ',i4,' WITH FORMAT: ',a)
    201 format(1x,///11x,a,' FOR LAYER',i4,/ &
               1X,'READING BINARY ON UNIT ',i4)
    202 format(1x,///11x,a,/ &
               1x,'READING BINARY ON UNIT ',i4)
    203 format(1x,///11x,a,' FOR CROSS SECTION',/ &
               1x,'READING BINARY ON UNIT ',i4)
    551 format(1x,i3,1x,60(1x,i1):/(5x,60(1x,i1)))
    552 format(1x,i3,1x,40(1x,i2):/(5x,40(1x,i2)))
    553 format(1x,i3,1x,30(1x,i3):/(5x,30(1x,i3)))
    554 format(1x,i3,1x,25(1x,i4):/(5x,25(1x,i4)))
    555 format(1x,i3,1x,20(1x,i5):/(5x,20(1x,i5)))
    556 format(1x,i3,1x,10(1x,i11):/(5x,10(1x,i11)))
    557 format(1x,i3,1x,25(1x,i2):/(5x,25(1x,i2)))
    558 format(1x,i3,1x,15(1x,i4):/(5x,15(1x,i4)))
    559 format(1x,i3,1x,10(1x,i6):/(5x,10(1x,i6)))
    601 format(' ERROR READING ARRAY CONTROL RECORD FOR ',a, &
               ' FOR LAYER',i4,':')
    602 format(' ERROR READING ARRAY CONTROL RECORD FOR ',a,':')
!     ------------------------------------------------------------------
!
!C1------READ ARRAY CONTROL RECORD AS CHARACTER DATA.
    call u9rdcom(in,iout,cntrl,ierr)
!
!C2------LOOK FOR ALPHABETIC WORD THAT INDICATES THAT THE RECORD IS FREE
!C2------FORMAT.  SET A FLAG SPECIFYING IF FREE FORMAT OR FIXED FORMAT.
    iclose=0
    ifree=1
    icol=1
    call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
    if (cntrl(istart:istop).eq.'CONSTANT') then
      locat=0
    elseif (cntrl(istart:istop).eq.'INTERNAL') then
      locat=in
    elseif (cntrl(istart:istop).eq.'EXTERNAL') then
      write(ermsg, *) 'Error. EXTERNAL arrays not supported in MODFLOW 6.'
      call store_error(ermsg)
      write(ermsg, *) 'Use OPEN/CLOSE instead.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
    elseif (cntrl(istart:istop).eq.'OPEN/CLOSE') then
      call urword(cntrl,icol,istart,istop,0,n,r,iout,in)
      fname=cntrl(istart:istop)
      locat=nunopn
      iclose=1
      ftype = 'OPEN/CLOSE'
    else
      write(ermsg, *) 'ERROR READING CONTROL RECORD FOR ' //                 &
                       trim(adjustl(aname))
      call store_error(ermsg)
      call store_error(trim(adjustl(cntrl)))
      write(ermsg, *) 'Use CONSTANT, INTERNAL, or OPEN/CLOSE.'
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
!
!C2A-----DID NOT FIND A RECOGNIZED WORD, SO NOT USING FREE FORMAT.
!C2A-----READ THE CONTROL RECORD THE ORIGINAL WAY.
      ifree=0
      read(cntrl,1,err=600) locat,iconst,fmtin,iprn
    endif
!
!C3------FOR FREE FORMAT CONTROL RECORD, READ REMAINING FIELDS.
    if (ifree.ne.0) then
      call urword(cntrl,icol,istart,istop,2,iconst,r,iout,in)
      if (locat.ne.0) then
        call urword(cntrl,icol,istart,istop,1,n,r,iout,in)
        fmtin=cntrl(istart:istop)
        if (iclose.ne.0) then
          if (fmtin.eq.'(BINARY)') then
            call openfile(locat, iout, fname, ftype, fmtarg_opt=FORM,    &
                          accarg_opt=ACCESS)
          else
                call openfile(locat, iout, fname, 'OPEN/CLOSE')
             endif
          endif
          if (locat.gt.0 .and. fmtin.eq.'(BINARY)') locat=-locat
          call urword(cntrl,icol,istart,istop,2,iprn,r,iout,in)
       endif
    endif
!
!C4------TEST LOCAT TO SEE HOW TO DEFINE ARRAY VALUES.
    if (locat.eq.0) then
!
!C4A-----LOCAT=0; SET ALL ARRAY VALUES EQUAL TO ICONST. RETURN.
      do i=1,ii
        do j=1,jj
          ia(j,i)=iconst
        enddo
      enddo
      if (k.gt.0 .and. iout>0) write(iout,82) aname,iconst,k
      if (k.le.0 .and. iout>0) write(iout,83) aname,iconst
      return
    elseif (locat.gt.0) then
!
!C4B-----LOCAT>0; READ FORMATTED RECORDS USING FORMAT FMTIN.
      if (iout>0) then
        if (k.gt.0) then
           write(iout,94) aname,k,locat,fmtin
           elseif (k.eq.0) then
              write(iout,95) aname,locat,fmtin
           else
              write(iout,96) aname,locat,fmtin
           endif
        endif
        do i=1,ii
          if (fmtin.eq.'(FREE)') then
             read(locat,*) (ia(j,i),j=1,jj)
          else
             read(locat,fmtin) (ia(j,i),j=1,jj)
          endif
        enddo
      else
!
!C4C-----LOCAT<0; READ UNFORMATTED RECORD CONTAINING ARRAY VALUES.
        locat=-locat
        if (iout>0) then
          if (k.gt.0) then
            write(iout,201) aname,k,locat
          elseif (k.eq.0) then
            write(iout,202) aname,locat
          else
            write(iout,203) aname,locat
          endif
        endif
        read(locat)
        read(locat) ia
      endif
!
!C5------IF ICONST NOT ZERO THEN MULTIPLY ARRAY VALUES BY ICONST.
      if (iclose.ne.0) close(unit=locat)
      if (iconst.ne.0) then
        do i=1,ii
          do j=1,jj
            ia(j,i)=ia(j,i)*iconst
          enddo
        enddo
      endif
!
!C6------IF PRINT CODE (IPRN) <0 THEN RETURN.
      if (iprn.lt.0) return
!
!C7------PRINT COLUMN NUMBERS AT TOP OF PAGE.
      if (iout>0) then
        if (iprn.gt.9 .or. iprn.eq.0) iprn=6
        select case (iprn)
        case (401)
          call ucolno(1,jj,4,60,2,iout)
        case (402)
          call ucolno(1,jj,4,40,3,iout)
        case (403)
          call ucolno(1,jj,4,30,4,iout)
        case (404)
          call ucolno(1,jj,4,25,5,iout)
        case (405)
          call ucolno(1,jj,4,20,6,iout)
        case (406)
          call ucolno(1,jj,4,10,12,iout)
        case (407)
          call ucolno(1,jj,4,25,3,iout)
        case (408)
          call ucolno(1,jj,4,15,5,iout)
        case (409)
          call ucolno(1,jj,4,10,7,iout)
        end select
!
!C8------PRINT EACH ROW IN THE ARRAY.
        do i=1,ii
          select case (iprn)
          case (501)
  !----------------FORMAT 60I1
            write(iout,551) i,(ia(j,i),j=1,jj)
          case (502)
  !----------------FORMAT 40I2
            write(iout,552) i,(ia(j,i),j=1,jj)
          case (503)
  !----------------FORMAT 30I3
            write(iout,553) i,(ia(j,i),j=1,jj)
          case (504)
  !----------------FORMAT 25I4
            write(iout,554) i,(ia(j,i),j=1,jj)
          case (505)
  !----------------FORMAT 20I5
            write(iout,555) i,(ia(j,i),j=1,jj)
          case (506)
  !----------------FORMAT 10I11
            write(iout,556) i,(ia(j,i),j=1,jj)
          case (507)
  !----------------FORMAT 25I2
            write(iout,557) i,(ia(j,i),j=1,jj)
          case (508)
  !----------------FORMAT 15I4
            write(iout,558) i,(ia(j,i),j=1,jj)
          case (509)
  !----------------FORMAT 10I6
            write(iout,559) i,(ia(j,i),j=1,jj)
          end select
        enddo
      endif
!
!C9------RETURN
      return
!
!C10-----CONTROL RECORD ERROR.
  600 if (k.gt.0) then
        write(ermsg,601) aname,k
      else
        write(ermsg,602) aname
      endif
      call store_error(ermsg)
      write(ermsg,'(1x,a)') cntrl
      call store_error(ermsg)
      call store_error_unit(in)
      call ustop()
  end subroutine u2dint

  subroutine u3dint_1darg(iin,iout,nlay,nrow,ncol,neq,ival,cval)
    ! dummy
    integer(I4B), intent(in) :: iin
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: neq
    integer(I4B), dimension(ncol*nrow*nlay), intent(inout) :: ival
    character (len=24), intent(in) :: cval
    !
    call u3dint(iin, iout, nlay, nrow, ncol, neq, ival, cval)
    !
    return
  end subroutine u3dint_1darg

  subroutine u3dint(iin,iout,nlay,nrow,ncol,neq,ival,cval)
! ******************************************************************************
! Read three-dimensional integer array, consisting of multiple 2d arrays with
! array headers.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    integer(I4B), intent(in) :: iin
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: neq
    integer(I4B), dimension(ncol,nrow,nlay), intent(inout) :: ival
    character (len=24), intent(in) :: cval
    !local
    integer(I4B) :: k
    !functions
! ------------------------------------------------------------------------------
    !code
    do k = 1, nlay
      call u2dint(ival(:,:,k), cval, nrow, ncol, k, iin, iout)
    enddo
    return
  end subroutine u3dint

  subroutine u3ddbl_1darg(iin,iout,nlay,nrow,ncol,neq,rval,cval)
    ! dummy
    integer(I4B), intent(in) :: iin
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: neq
    real(DP), dimension(ncol*nrow*nlay), intent(inout) :: rval
    character (len=24), intent(in) :: cval
    !
    call u3ddbl(iin, iout, nlay, nrow, ncol, neq, rval, cval)
    !
    return
  end subroutine u3ddbl_1darg

  subroutine u3ddbl(iin,iout,nlay,nrow,ncol,neq,rval,cval)
! ******************************************************************************
! Read three-dimensional double-precision array, consisting of multiple 2d
! arrays with array headers.
! ******************************************************************************
!
!    SPECIFICATIONS:
! ------------------------------------------------------------------------------
    implicit none
    integer(I4B), intent(in) :: iin
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: neq
    real(DP), dimension(ncol,nrow,nlay), intent(inout) :: rval
    character (len=24), intent(in) :: cval
    !local
    integer(I4B) :: k
    !functions
! ------------------------------------------------------------------------------
    !code
!    istart = 1
!    istop = ncol * nrow
    do k = 1, nlay
      call u2ddbl(rval(:,:,k),cval,nrow,ncol,k,iin,iout)
!      istart = istop + 1
!      istop = istop + ncol * nrow
    enddo
    return
  end subroutine u3ddbl

end module ArrayReadersMF5Module
