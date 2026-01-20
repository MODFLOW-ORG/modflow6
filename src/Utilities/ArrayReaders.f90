module ArrayReadersModule

  use ConstantsModule, only: DONE, LINELENGTH, LENBIGLINE, LENBOUNDNAME, &
                             NAMEDBOUNDFLAG, LINELENGTH, DZERO, MAXCHARLEN, &
                             DZERO
  use InputOutputModule, only: openfile, u9rdcom, urword, ucolno, ulaprw, &
                               BuildFixedFormat, BuildFloatFormat, &
                               BuildIntFormat
  use KindModule, only: DP, I4B, LGP
  use OpenSpecModule, only: ACCESS, FORM
  use SimModule, only: store_error, store_error_unit
  use SimVariablesModule, only: errmsg

  implicit none

  private
  public :: ReadArray
  public :: read_binary_header
  public :: check_binary_filesize
  public :: BINARY_INT_BYTES
  public :: BINARY_DOUBLE_BYTES
  public :: BINARY_HEADER_BYTES

  integer(I4B), parameter :: BINARY_CHAR_BYTES = 1
  integer(I4B), parameter :: BINARY_INT_BYTES = 4
  integer(I4B), parameter :: BINARY_DOUBLE_BYTES = 8
  integer(I4B), parameter :: BINARY_STRLEN = 16
  integer(I4B), parameter :: BINARY_HEADER_BYTES = &
                             (5 * BINARY_INT_BYTES) + & !< kstp, kper, msize1, msize2, msize3
                             (2 * BINARY_DOUBLE_BYTES) + & !< pertim, totim
                             (BINARY_STRLEN * BINARY_CHAR_BYTES) !< array text

  interface ReadArray
    module procedure &
      read_array_int1d, &
      read_array_int2d, &
      read_array_int3d, &
      read_array_dbl1d, &
      read_array_dbl2d, &
      read_array_dbl3d, &
      read_array_dbl1d_layered, &
      read_array_int1d_layered, &
      read_array_dbl3d_all, &
      read_array_int3d_all
  end interface ReadArray

  ! Integer readers
  ! read_array_int1d(iu, iarr, aname, ndim, jj, iout, k)
  ! read_array_int1d_layered(iu, iarr, aname, ndim, ncol, nrow, nlay, nval, iout, k1, k2)
  ! read_array_int2d(iu, iarr, aname, ndim, jj, ii, iout, k)
  ! read_array_int3d(iu, iarr, aname, ndim, ncol, nrow, nlay, iout, k1, k2)
  ! read_array_int3d_all(iu, iarr, aname, ndim, nvals, iout)
  !
  ! Floating-point readers
  ! read_array_dbl1d(iu, darr, aname, ndim, jj, iout, k)
  ! read_array_dbl1d_layered(iu, darr, aname, ndim, ncol, nrow, nlay, nval, iout, k1, k2)
  ! read_array_dbl2d(iu, darr, aname, ndim, jj, ii, iout, k)
  ! read_array_dbl3d(iu, darr, aname, ndim, ncol, nrow, nlay, iout, k1, k2)
  ! read_array_dbl3d_all(iu, darr, aname, ndim, nvals, iout)

contains

  ! -- Procedures that are part of ReadArray interface (integer data)

  subroutine read_array_int1d(iu, iarr, aname, ndim, jj, iout, k)
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout
    integer(I4B), intent(in) :: jj
    integer(I4B), dimension(jj), intent(inout) :: iarr
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: ndim ! dis%ndim
    integer(I4B), intent(in) :: k ! layer number; 0 to not print
    ! -- local
    logical(LGP) :: isok
    integer(I4B) :: iclose, iconst, iprn, j, locat, ncpl, ndig
    integer(I4B) :: nval, nvalt
    logical :: prowcolnum
    character(len=100) :: prfmt
    integer(I4B) :: istat
    character(len=30) :: arrname
    character(len=MAXCHARLEN) :: ermsgr
    ! -- formats
2   format(/, 1x, a, ' = ', i0, ' FOR LAYER ', i0)
3   format(/, 1x, a, ' = ', i0)
    !
    ! -- Read array control record.
    call read_control_int(iu, iout, aname, locat, iconst, iclose, iprn)
    !
    ! -- Read or assign array data.
    if (locat == 0) then
      ! -- Assign constant
      do j = 1, jj
        iarr(j) = iconst
      end do
      if (iout > 0) then
        if (k > 0) then
          write (iout, 2) trim(aname), iconst, k
        else
          write (iout, 3) trim(aname), iconst
        end if
      end if
    elseif (locat > 0) then
      ! -- Read data as text
      read (locat, *, iostat=istat, iomsg=ermsgr) (iarr(j), j=1, jj)
      if (istat /= 0) then
        arrname = adjustl(aname)
        errmsg = "Error reading data for array '"//trim(arrname)// &
                 "'. "//trim(adjustl(ermsgr))
        call store_error(errmsg)
        call store_error_unit(locat)
      end if
      do j = 1, jj
        iarr(j) = iarr(j) * iconst
      end do
      if (iclose == 1) then
        close (locat)
      end if
    else
      ! -- Read data as binary
      locat = -locat
      nvalt = 0
      do
        call read_binary_header(locat, iout, aname, nval)
        isok = check_binary_size(nval, nvalt, size(iarr), aname, locat)
        if (isok .EQV. .FALSE.) exit
        read (locat, iostat=istat, iomsg=ermsgr) &
          (iarr(j), j=nvalt + 1, nvalt + nval)
        if (istat /= 0) then
          arrname = adjustl(aname)
          errmsg = "Error reading data for array '"//trim(arrname)// &
                   "'. "//trim(adjustl(ermsgr))
          call store_error(errmsg)
          call store_error_unit(locat)
        end if
        nvalt = nvalt + nval
        if (nvalt == size(iarr)) exit
      end do
      !
      ! -- multiply array by constant
      do j = 1, jj
        iarr(j) = iarr(j) * iconst
      end do
      !
      ! -- close the file
      if (iclose == 1) then
        close (locat)
      end if
    end if
    !
    ! -- Print array if requested.
    if (iprn >= 0 .and. locat /= 0) then
      prowcolnum = (ndim == 3)
      call build_format_int(iprn, prfmt, prowcolnum, ncpl, ndig)
      call print_array_int(iarr, aname, iout, jj, 1, k, prfmt, ncpl, ndig, &
                           prowcolnum)
    end if
  end subroutine read_array_int1d

  subroutine read_array_int2d(iu, iarr, aname, ndim, jj, ii, iout, k)
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout
    integer(I4B), intent(in) :: jj, ii
    integer(I4B), dimension(jj, ii), intent(inout) :: iarr
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: ndim ! dis%ndim
    integer(I4B), intent(in) :: k ! layer number; 0 to not print
    ! -- local
    logical(LGP) :: isok
    integer(I4B) :: i, iclose, iconst, iprn, j, locat, ncpl, ndig
    integer(I4B) :: nval
    logical :: prowcolnum
    character(len=100) :: prfmt
    integer(I4B) :: istat
    character(len=30) :: arrname
    character(len=MAXCHARLEN) :: ermsgr
    ! -- formats
2   format(/, 1x, a, ' = ', i0, ' FOR LAYER ', i0)
3   format(/, 1x, a, ' = ', i0)
    !
    ! -- Read array control record.
    call read_control_int(iu, iout, aname, locat, iconst, iclose, iprn)
    !
    ! -- Read or assign array data.
    if (locat == 0) then
      ! -- Assign constant
      do i = 1, ii
        do j = 1, jj
          iarr(j, i) = iconst
        end do
      end do
      if (iout > 0) then
        if (k > 0) then
          write (iout, 2) trim(aname), iconst, k
        else
          write (iout, 3) trim(aname), iconst
        end if
      end if
    elseif (locat > 0) then
      ! -- Read data as text
      do i = 1, ii
        read (locat, *, iostat=istat, iomsg=ermsgr) (iarr(j, i), j=1, jj)
        if (istat /= 0) then
          arrname = adjustl(aname)
          errmsg = "Error reading data for array '"//trim(arrname)// &
                   "'. "//trim(adjustl(ermsgr))
          call store_error(errmsg)
          call store_error_unit(locat)
        end if
        do j = 1, jj
          iarr(j, i) = iarr(j, i) * iconst
        end do
      end do
      if (iclose == 1) then
        close (locat)
      end if
    else
      ! -- Read data as binary
      locat = -locat
      call read_binary_header(locat, iout, aname, nval)
      isok = check_binary_size(nval, 0, size(iarr), aname, locat)
      if (isok) then
        do i = 1, ii
          read (locat, iostat=istat, iomsg=ermsgr) (iarr(j, i), j=1, jj)
          if (istat /= 0) then
            arrname = adjustl(aname)
            errmsg = "Error reading data for array '"//trim(arrname)// &
                     "'. "//trim(adjustl(ermsgr))
            call store_error(errmsg)
            call store_error_unit(locat)
          end if
          do j = 1, jj
            iarr(j, i) = iarr(j, i) * iconst
          end do
        end do
      end if
      if (iclose == 1) then
        close (locat)
      end if
    end if
    !
    ! -- Print array if requested.
    if (iprn >= 0 .and. locat /= 0) then
      prowcolnum = (ndim == 3)
      call build_format_int(iprn, prfmt, prowcolnum, ncpl, ndig)
      call print_array_int(iarr, aname, iout, jj, ii, k, prfmt, ncpl, &
                           ndig, prowcolnum)
    end if
  end subroutine read_array_int2d

  subroutine read_array_int3d(iu, iarr, aname, ndim, ncol, nrow, nlay, iout, &
                              k1, k2)
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: ndim
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: k1, k2
    integer(I4B), dimension(ncol, nrow, nlay), intent(inout) :: iarr
    character(len=*), intent(in) :: aname
    ! -- local
    integer(I4B) :: k, kk
    do k = k1, k2
      if (k <= 0) then
        kk = 1
      else
        kk = k
      end if
      call read_array_int2d(iu, iarr(:, :, kk), aname, ndim, ncol, nrow, iout, k)
    end do
  end subroutine read_array_int3d

  subroutine read_array_int3d_all(iu, iarr, aname, ndim, nvals, iout)
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: ndim
    integer(I4B), intent(in) :: nvals
    integer(I4B), dimension(nvals, 1, 1), intent(inout) :: iarr
    character(len=*), intent(in) :: aname
    ! -- local
    !
    call read_array_int1d(iu, iarr, aname, ndim, nvals, iout, 0)
  end subroutine read_array_int3d_all

  subroutine read_array_int1d_layered(iu, iarr, aname, ndim, ncol, nrow, &
                                      nlay, nval, iout, k1, k2)
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout
    integer(I4B), intent(in) :: ncol, nrow, nlay, nval
    integer(I4B), dimension(nval), intent(inout) :: iarr
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: ndim ! dis%ndim
    integer(I4B), intent(in) :: k1, k2
    ! -- local
    !
    call read_array_int3d(iu, iarr, aname, ndim, ncol, nrow, nlay, iout, k1, k2)
  end subroutine read_array_int1d_layered

  ! -- Procedures that are part of ReadArray interface (floating-point data)

  subroutine read_array_dbl1d(iu, darr, aname, ndim, jj, iout, k)
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout
    integer(I4B), intent(in) :: jj
    real(DP), dimension(jj), intent(inout) :: darr
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: ndim ! dis%ndim
    integer(I4B), intent(in) :: k ! layer number; 0 to not print
    ! -- local
    logical(LGP) :: isok
    integer(I4B) :: j, iclose, iprn, locat, ncpl, ndig
    real(DP) :: cnstnt
    logical :: prowcolnum
    character(len=100) :: prfmt
    integer(I4B) :: istat
    integer(I4B) :: nvalt, nval
    character(len=30) :: arrname
    character(len=MAXCHARLEN) :: ermsgr
    ! -- formats
2   format(/, 1x, a, ' = ', g14.7, ' FOR LAYER ', i0)
3   format(/, 1x, a, ' = ', g14.7)
    !
    ! -- Read array control record.
    call read_control_dbl(iu, iout, aname, locat, cnstnt, iclose, iprn)
    !
    ! -- Read or assign array data.
    if (locat == 0) then
      ! -- Assign constant
      do j = 1, jj
        darr(j) = cnstnt
      end do
      if (iout > 0) then
        if (k > 0) then
          write (iout, 2) trim(aname), cnstnt, k
        else
          write (iout, 3) trim(aname), cnstnt
        end if
      end if
    elseif (locat > 0) then
      ! -- Read data as text
      read (locat, *, iostat=istat, iomsg=ermsgr) (darr(j), j=1, jj)
      if (istat /= 0) then
        arrname = adjustl(aname)
        errmsg = "Error reading data for array '"// &
                 trim(adjustl(arrname))//"'. "//trim(adjustl(ermsgr))
        call store_error(errmsg)
        call store_error_unit(locat)
      end if
      do j = 1, jj
        darr(j) = darr(j) * cnstnt
      end do
      if (iclose == 1) then
        close (locat)
      end if
    else
      ! -- Read data as binary
      locat = -locat
      nvalt = 0
      do
        call read_binary_header(locat, iout, aname, nval)
        isok = check_binary_size(nval, nvalt, size(darr), aname, locat)
        if (isok .EQV. .FALSE.) exit
        read (locat, iostat=istat, iomsg=ermsgr) &
          (darr(j), j=nvalt + 1, nvalt + nval)
        if (istat /= 0) then
          arrname = adjustl(aname)
          errmsg = "Error reading data for array '"// &
                   trim(adjustl(arrname))//"'. "//trim(adjustl(ermsgr))
          call store_error(errmsg)
          call store_error_unit(locat)
        end if
        nvalt = nvalt + nval
        if (nvalt == size(darr)) exit
      end do
      !
      ! -- multiply entire array by constant
      do j = 1, jj
        darr(j) = darr(j) * cnstnt
      end do
      !
      ! -- close the file
      if (iclose == 1) then
        close (locat)
      end if
    end if
    !
    ! -- Print array if requested.
    if (iprn >= 0 .and. locat /= 0) then
      prowcolnum = (ndim == 3)
      call build_format_dbl(iprn, prfmt, prowcolnum, ncpl, ndig)
      call print_array_dbl(darr, aname, iout, jj, 1, k, prfmt, ncpl, ndig, &
                           prowcolnum)
    end if
  end subroutine read_array_dbl1d

  subroutine read_array_dbl2d(iu, darr, aname, ndim, jj, ii, iout, k)
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout
    integer(I4B), intent(in) :: jj, ii
    real(DP), dimension(jj, ii), intent(inout) :: darr
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: ndim ! dis%ndim
    integer(I4B), intent(in) :: k ! layer number; 0 to not print
    ! -- local
    logical(LGP) :: isok
    integer(I4B) :: i, iclose, iprn, j, locat, ncpl, ndig
    integer(I4B) :: nval
    real(DP) :: cnstnt
    logical :: prowcolnum
    character(len=100) :: prfmt
    integer(I4B) :: istat
    character(len=30) :: arrname
    character(len=MAXCHARLEN) :: ermsgr
    ! -- formats
2   format(/, 1x, a, ' = ', g14.7, ' FOR LAYER ', i0)
3   format(/, 1x, a, ' = ', g14.7)
    !
    ! -- Read array control record.
    call read_control_dbl(iu, iout, aname, locat, cnstnt, iclose, iprn)
    !
    ! -- Read or assign array data.
    if (locat == 0) then
      ! -- Assign constant
      do i = 1, ii
        do j = 1, jj
          darr(j, i) = cnstnt
        end do
      end do
      if (iout > 0) then
        if (k > 0) then
          write (iout, 2) trim(aname), cnstnt, k
        else
          write (iout, 3) trim(aname), cnstnt
        end if
      end if
    elseif (locat > 0) then
      ! -- Read data as text
      do i = 1, ii
        read (locat, *, iostat=istat, iomsg=ermsgr) (darr(j, i), j=1, jj)
        if (istat /= 0) then
          arrname = adjustl(aname)
          errmsg = "Error reading data for array '"// &
                   trim(adjustl(arrname))//"'. "//trim(adjustl(ermsgr))
          call store_error(errmsg)
          call store_error_unit(locat)
        end if
        do j = 1, jj
          darr(j, i) = darr(j, i) * cnstnt
        end do
      end do
      if (iclose == 1) then
        close (locat)
      end if
    else
      ! -- Read data as binary
      locat = -locat
      call read_binary_header(locat, iout, aname, nval)
      isok = check_binary_size(nval, 0, size(darr), aname, locat)
      if (isok) then
        do i = 1, ii
          read (locat, iostat=istat, iomsg=ermsgr) (darr(j, i), j=1, jj)
          if (istat /= 0) then
            arrname = adjustl(aname)
            errmsg = "Error reading data for array '"// &
                     trim(adjustl(arrname))//"'. "//trim(adjustl(ermsgr))
            call store_error(errmsg)
            call store_error_unit(locat)
          end if
          do j = 1, jj
            darr(j, i) = darr(j, i) * cnstnt
          end do
        end do
      end if
      if (iclose == 1) then
        close (locat)
      end if
    end if
    !
    ! -- Print array if requested.
    if (iprn >= 0 .and. locat /= 0) then
      prowcolnum = (ndim == 3)
      call build_format_dbl(iprn, prfmt, prowcolnum, ncpl, ndig)
      call print_array_dbl(darr, aname, iout, jj, ii, k, prfmt, ncpl, &
                           ndig, prowcolnum)
    end if
  end subroutine read_array_dbl2d

  subroutine read_array_dbl3d(iu, darr, aname, ndim, ncol, nrow, nlay, iout, &
                              k1, k2)
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: ndim
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(in) :: k1, k2
    real(DP), dimension(ncol, nrow, nlay), intent(inout) :: darr
    character(len=*), intent(in) :: aname
    ! -- local
    integer(I4B) :: k, kk
    !
    do k = k1, k2
      if (k <= 0) then
        kk = 1
      else
        kk = k
      end if
      call read_array_dbl2d(iu, darr(:, :, kk), aname, ndim, ncol, nrow, iout, k)
    end do
  end subroutine read_array_dbl3d

  subroutine read_array_dbl3d_all(iu, darr, aname, ndim, nvals, iout)
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    integer(I4B), intent(in) :: ndim
    integer(I4B), intent(in) :: nvals
    real(DP), dimension(nvals, 1, 1), intent(inout) :: darr
    character(len=*), intent(in) :: aname
    ! -- local
    !
    call read_array_dbl1d(iu, darr, aname, ndim, nvals, iout, 0)
  end subroutine read_array_dbl3d_all

  subroutine read_array_dbl1d_layered(iu, darr, aname, ndim, ncol, nrow, &
                                      nlay, nval, iout, k1, k2)
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout
    integer(I4B), intent(in) :: ncol, nrow, nlay, nval
    real(DP), dimension(nval), intent(inout) :: darr
    character(len=*), intent(in) :: aname
    integer(I4B), intent(in) :: ndim ! dis%ndim
    integer(I4B), intent(in) :: k1, k2
    ! -- local
    !
    call read_array_dbl3d(iu, darr, aname, ndim, ncol, nrow, nlay, iout, k1, k2)
  end subroutine read_array_dbl1d_layered

  ! -- Utility procedures

  subroutine read_control_int(iu, iout, aname, locat, iconst, &
                              iclose, iprn)
    ! Read an array-control record for an integer array.
    ! Open an input file if needed.
    ! If CONSTANT is specified in input, locat is returned as 0.
    ! If (BINARY) is specified, locat is returned as the negative of
    ! the unit number opened for binary read.
    ! If OPEN/CLOSE is specified, iclose is returned as 1, otherwise 0.
    ! -- dummy
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: aname
    integer(I4B), intent(out) :: locat
    integer(I4B), intent(out) :: iconst
    integer(I4B), intent(out) :: iclose
    integer(I4B), intent(out) :: iprn
    ! -- local
    integer(I4B) :: icol, icol1, istart, istop, n
    real(DP) :: r
    character(len=MAXCHARLEN) :: fname
    character(len=:), allocatable :: line
    !
    ! -- Read CONSTANT, INTERNAL, or OPEN/CLOSE from array control record.
    call read_control_1(iu, iout, aname, locat, iclose, line, icol, fname)
    if (locat == 0) then
      ! CONSTANT was found -- read value and return
      call urword(line, icol, istart, istop, 2, iconst, r, iout, iu)
      iprn = -1
      return
    end if
    icol1 = icol
    iconst = 1
    !
    ! -- Read FACTOR option from array control record.
    call urword(line, icol, istart, istop, 1, n, r, iout, iu)
    if (line(istart:istop) == 'FACTOR') then
      call urword(line, icol, istart, istop, 2, iconst, r, iout, iu)
      if (iconst == 0) iconst = 1
    else
      icol = icol1
    end if
    !
    ! -- Read (BINARY) and IPRN options from array control record,
    !    and open an OPEN/CLOSE file if specified.
    call read_control_2(iu, iout, fname, line, icol, locat, iclose, iprn)
  end subroutine read_control_int

  subroutine read_control_dbl(iu, iout, aname, locat, cnstnt, &
                              iclose, iprn)
    ! Read an array-control record for a double-precision array.
    ! Open an input file if needed.
    ! If CONSTANT is specified in input, locat is returned as 0.
    ! If (BINARY) is specified, locat is returned as the negative of
    ! the unit number opened for binary read.
    ! If OPEN/CLOSE is specified, iclose is returned as 1, otherwise 0.
    ! -- dummy
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: aname
    integer(I4B), intent(out) :: locat
    real(DP), intent(out) :: cnstnt
    integer(I4B), intent(out) :: iclose
    integer(I4B), intent(out) :: iprn
    !
    ! -- local
    integer(I4B) :: icol, icol1, istart, istop, n
    real(DP) :: r
    character(len=MAXCHARLEN) :: fname
    character(len=:), allocatable :: line
    !
    ! -- Read CONSTANT, INTERNAL, or OPEN/CLOSE from array control record.
    call read_control_1(iu, iout, aname, locat, iclose, line, icol, fname)
    if (locat == 0) then
      ! CONSTANT was found -- read value and return
      call urword(line, icol, istart, istop, 3, n, cnstnt, iout, iu)
      iprn = -1
      return
    end if
    icol1 = icol
    cnstnt = DONE
    !
    ! -- Read FACTOR option from array control record.
    call urword(line, icol, istart, istop, 1, n, r, iout, iu)
    if (line(istart:istop) == 'FACTOR') then
      call urword(line, icol, istart, istop, 3, n, cnstnt, iout, iu)
      if (cnstnt == DZERO) cnstnt = DONE
    else
      icol = icol1
    end if
    !
    ! -- Read (BINARY) and IPRN options from array control record,
    !    and open an OPEN/CLOSE file if specified.
    call read_control_2(iu, iout, fname, line, icol, locat, iclose, iprn)
  end subroutine read_control_dbl

  subroutine read_control_1(iu, iout, aname, locat, iclose, line, icol, fname)
    use SimModule, only: ustop
    ! -- Read CONSTANT, INTERNAL, or OPEN/CLOSE from array control record.
    ! -- dummy
    integer(I4B), intent(in) :: iu
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: aname
    integer(I4B), intent(out) :: locat
    integer(I4B), intent(out) :: iclose
    character(len=:), allocatable, intent(inout) :: line
    integer(I4B), intent(inout) :: icol
    character(len=*), intent(inout) :: fname

    ! -- local
    integer(I4B) :: istart, istop, n
    integer(I4B) :: ierr
    real(DP) :: r
    !
    ! -- Read array control record.  Any future refactoring
    !    should use the LongLineReader here instead of u9rdcom
    call u9rdcom(iu, iout, line, ierr)
    !
    iclose = 0
    icol = 1
    ! -- Read first token of array control record.
    call urword(line, icol, istart, istop, 1, n, r, iout, iu)
    if (line(istart:istop) .eq. 'CONSTANT') then
      locat = 0
    elseif (line(istart:istop) .eq. 'INTERNAL') then
      locat = iu
    elseif (line(istart:istop) .eq. 'OPEN/CLOSE') then
      call urword(line, icol, istart, istop, 0, n, r, iout, iu)
      fname = line(istart:istop)
      locat = -1
      iclose = 1
    else
      errmsg = 'READING CONTROL RECORD FOR '// &
               trim(adjustl(aname))//"'. "// &
               'Use CONSTANT, INTERNAL, or OPEN/CLOSE.'
      call store_error(errmsg)
      call store_error_unit(iu)
    end if
  end subroutine read_control_1

  subroutine read_control_2(iu, iout, fname, line, icol, &
                            locat, iclose, iprn)
    ! -- Read (BINARY) and IPRN options from array control record,
    !    and open an OPEN/CLOSE file if specified.
    ! -- dummy
    integer(I4B), intent(in) :: iu, iout, iclose
    character(len=*), intent(in) :: fname
    character(len=*), intent(inout) :: line
    integer(I4B), intent(inout) :: icol, iprn, locat
    ! -- local
    integer(I4B) :: i, n, istart, istop, lenkey
    real(DP) :: r
    character(len=MAXCHARLEN) :: keyword
    logical :: binary
    !
    iprn = -1 ! Printing is turned off by default
    binary = .false.
    !
    if (locat .ne. 0) then
      ! -- CONSTANT has not been specified; array data will be read.
      ! -- Read at most two options.
      do i = 1, 2
        call urword(line, icol, istart, istop, 1, n, r, iout, iu)
        keyword = line(istart:istop)
        lenkey = len_trim(keyword)
        select case (keyword)
        case ('(BINARY)')
          if (iclose == 0) then
            errmsg = '"(BINARY)" option for array input is valid only if'// &
                     ' OPEN/CLOSE is also specified.'
            call store_error(errmsg)
            call store_error_unit(iu)
          end if
          binary = .true.
        case ('IPRN')
          ! -- Read IPRN value
          call urword(line, icol, istart, istop, 2, iprn, r, iout, iu)
          exit
        case ('')
          exit
        case default
          errmsg = 'Invalid option found in array-control record: "' &
                   //trim(keyword)//'"'
          call store_error(errmsg)
          call store_error_unit(iu)
        end select
      end do
      !
      if (iclose == 0) then
        ! -- Array data will be read from current input file.
        locat = iu
      else
        ! -- Open the OPEN\CLOSE file
        if (binary) then
          call openfile(locat, iout, fname, 'OPEN/CLOSE', fmtarg_opt=FORM, &
                        accarg_opt=ACCESS)
          locat = -locat
        else
          call openfile(locat, iout, fname, 'OPEN/CLOSE')
        end if
      end if
    end if
  end subroutine read_control_2

  subroutine build_format_int(iprn, prfmt, prowcolnum, ncpl, ndig)
    ! -- Build a print format for integers based on IPRN.
    ! -- dummy
    integer(I4B), intent(inout) :: iprn
    character(len=*), intent(out) :: prfmt
    logical, intent(in) :: prowcolnum
    integer(I4B), intent(out) :: ncpl, ndig
    ! -- local
    integer(I4B) :: nwidp
    !
    if (iprn < 0) then
      prfmt = ''
      return
    end if
    !
    if (iprn > 9) iprn = 0
    !
    select case (iprn)
    case (0)
      ncpl = 10
      nwidp = 11
    case (1)
      ncpl = 60
      nwidp = 1
    case (2)
      ncpl = 40
      nwidp = 2
    case (3)
      ncpl = 30
      nwidp = 3
    case (4)
      ncpl = 25
      nwidp = 4
    case (5)
      ncpl = 20
      nwidp = 5
    case (6)
      ncpl = 10
      nwidp = 11
    case (7)
      ncpl = 25
      nwidp = 2
    case (8)
      ncpl = 15
      nwidp = 4
    case (9)
      ncpl = 19
      nwidp = 6
    end select
    !
    call BuildIntFormat(ncpl, nwidp, prfmt, prowcolnum)
    ndig = nwidp + 1
  end subroutine build_format_int

  subroutine build_format_dbl(iprn, prfmt, prowcolnum, ncpl, ndig)
    ! -- Build a print format for reals based on IPRN.
    ! -- dummy
    integer(I4B), intent(inout) :: iprn
    character(len=*), intent(out) :: prfmt
    logical, intent(in) :: prowcolnum
    integer(I4B), intent(out) :: ncpl, ndig
    ! -- local
    integer(I4B) :: nwidp
    character(len=1) :: editdesc
    !
    if (iprn < 0) then
      prfmt = ''
      return
    end if
    !
    if (iprn > 21) iprn = 0
    !
    select case (iprn)
    case (0)
      ncpl = 10
      editdesc = 'G'
      nwidp = 11
      ndig = 4
    case (1)
      ncpl = 11
      editdesc = 'G'
      nwidp = 10
      ndig = 3
    case (2)
      ncpl = 9
      editdesc = 'G'
      nwidp = 13
      ndig = 6
    case (3)
      ncpl = 15
      editdesc = 'F'
      nwidp = 7
      ndig = 1
    case (4)
      ncpl = 15
      editdesc = 'F'
      nwidp = 7
      ndig = 2
    case (5)
      ncpl = 15
      editdesc = 'F'
      nwidp = 7
      ndig = 3
    case (6)
      ncpl = 15
      editdesc = 'F'
      nwidp = 7
      ndig = 4
    case (7)
      ncpl = 20
      editdesc = 'F'
      nwidp = 5
      ndig = 0
    case (8)
      ncpl = 20
      editdesc = 'F'
      nwidp = 5
      ndig = 1
    case (9)
      ncpl = 20
      editdesc = 'F'
      nwidp = 5
      ndig = 2
    case (10)
      ncpl = 20
      editdesc = 'F'
      nwidp = 5
      ndig = 3
    case (11)
      ncpl = 20
      editdesc = 'F'
      nwidp = 5
      ndig = 4
    case (12)
      ncpl = 10
      editdesc = 'G'
      nwidp = 11
      ndig = 4
    case (13)
      ncpl = 10
      editdesc = 'F'
      nwidp = 6
      ndig = 0
    case (14)
      ncpl = 10
      editdesc = 'F'
      nwidp = 6
      ndig = 1
    case (15)
      ncpl = 10
      editdesc = 'F'
      nwidp = 6
      ndig = 2
    case (16)
      ncpl = 10
      editdesc = 'F'
      nwidp = 6
      ndig = 3
    case (17)
      ncpl = 10
      editdesc = 'F'
      nwidp = 6
      ndig = 4
    case (18)
      ncpl = 10
      editdesc = 'F'
      nwidp = 6
      ndig = 5
    case (19)
      ncpl = 5
      editdesc = 'G'
      nwidp = 12
      ndig = 5
    case (20)
      ncpl = 6
      editdesc = 'G'
      nwidp = 11
      ndig = 4
    case (21)
      ncpl = 7
      editdesc = 'G'
      nwidp = 9
      ndig = 2
    end select
    !
    if (editdesc == 'F') then
      call BuildFixedFormat(ncpl, nwidp, ndig, prfmt, prowcolnum)
    else
      call BuildFloatFormat(ncpl, nwidp, ndig, editdesc, prfmt, prowcolnum)
    end if
    !
    ndig = nwidp + 1
  end subroutine build_format_dbl

  subroutine print_array_int(iarr, aname, iout, jj, ii, k, prfmt, &
                             ncpl, ndig, prowcolnum)
    ! -- dummy
    integer(I4B), intent(in) :: iout, jj, ii, k
    integer(I4B), intent(in) :: ncpl ! # values to print per line
    integer(I4B), intent(in) :: ndig ! # characters in each field
    integer(I4B), dimension(jj, ii), intent(in) :: iarr ! Integer array to be printed
    character(len=*), intent(in) :: aname ! Array name
    character(len=*), intent(in) :: prfmt ! Print format, no row #
    logical, intent(in) :: prowcolnum ! Print row & column numbers
    ! -- local
    integer(I4B) :: i, j
    ! -- formats
2   format(/, 1x, a, 1x, 'FOR LAYER ', i0)
3   format(/, 1x, a)
    !
    if (iout <= 0) return
    !
    ! -- Write name of array
    if (k > 0) then
      write (iout, 2) trim(aname), k
    else
      write (iout, 3) trim(aname)
    end if
    !
    ! -- Write array
    if (prowcolnum) then
      ! -- Write column/node numbers
      call ucolno(1, jj, 4, ncpl, ndig, iout)
      !
      ! -- Write array values, including row numbers
      do i = 1, ii
        write (iout, prfmt) i, (iarr(j, i), j=1, jj)
      end do
    else
      if (ii > 1) then
        errmsg = 'Program error printing array '//trim(aname)// &
                 ': ii > 1 when prowcolnum is false.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
      !
      ! -- Write array values, without row numbers
      write (iout, prfmt) (iarr(j, 1), j=1, jj)
    end if
  end subroutine print_array_int

  subroutine print_array_dbl(darr, aname, iout, jj, ii, k, prfmt, &
                             ncpl, ndig, prowcolnum)
    ! -- dummy
    integer(I4B), intent(in) :: iout, jj, ii, k
    integer(I4B), intent(in) :: ncpl ! # values to print per line
    integer(I4B), intent(in) :: ndig ! # characters in each field
    real(DP), dimension(jj, ii), intent(in) :: darr ! Real array to be printed
    character(len=*), intent(in) :: aname ! Array name
    character(len=*), intent(in) :: prfmt ! Print format, no row #
    logical, intent(in) :: prowcolnum ! Print row & column numbers
    ! -- local
    integer(I4B) :: i, j
    ! -- formats
2   format(/, 1x, a, 1x, 'FOR LAYER ', i0)
3   format(/, 1x, a)
    !
    if (iout <= 0) return
    !
    ! -- Write name of array
    if (k > 0) then
      write (iout, 2) trim(aname), k
    else
      write (iout, 3) trim(aname)
    end if
    !
    ! -- Write array
    if (prowcolnum) then
      ! -- Write column/node numbers
      call ucolno(1, jj, 4, ncpl, ndig, iout)
      !
      ! -- Write array values, including row numbers
      do i = 1, ii
        write (iout, prfmt) i, (darr(j, i), j=1, jj)
      end do
    else
      if (ii > 1) then
        errmsg = 'Program error printing array '//trim(aname)// &
                 ': ii > 1 when prowcolnum is false.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
      !
      ! -- Write array values, without row numbers
      write (iout, prfmt) (darr(j, 1), j=1, jj)
    end if
  end subroutine print_array_dbl

  subroutine read_binary_header(locat, iout, arrname, nval)
    ! -- dummy
    integer(I4B), intent(in) :: locat
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: arrname
    integer, intent(out) :: nval
    ! -- local
    integer(I4B) :: istat
    integer(I4B) :: kstp, kper, m1, m2, m3
    real(DP) :: pertim, totim
    character(len=BINARY_STRLEN) :: text
    character(len=MAXCHARLEN) :: ermsgr
    character(len=*), parameter :: fmthdr = &
      "(/,1X,'HEADER FROM BINARY FILE HAS FOLLOWING ENTRIES',&
       &/,4X,'KSTP: ',I0,'  KPER: ',I0,&
       &/,4x,'PERTIM: ',G0,'  TOTIM: ',G0,&
       &/,4X,'TEXT: ',A,&
       &/,4X,'MSIZE 1: ',I0,'  MSIZE 2: ',I0,'  MSIZE 3: ',I0)"
    !
    ! -- Read the header line from the binary file
    read (locat, iostat=istat, iomsg=ermsgr) kstp, kper, pertim, totim, text, &
      m1, m2, m3
    !
    ! -- Check for errors
    if (istat /= 0) then
      errmsg = "Error reading data for array '"//adjustl(trim(arrname))// &
               "'. "//trim(adjustl(ermsgr))
      call store_error(errmsg)
      call store_error_unit(locat)
    end if
    !
    ! -- Write message about the binary header
    if (iout > 0) then
      write (iout, fmthdr) kstp, kper, pertim, totim, text, m1, m2, m3
    end if
    !
    ! -- Assign the number of values that follow the header
    nval = m1 * m2
  end subroutine read_binary_header

  subroutine check_binary_filesize(locat, expected_size, arrname)
    ! -- dummy
    integer(I4B), intent(in) :: locat
    integer(I4B), intent(in) :: expected_size
    character(len=*), intent(in) :: arrname
    ! -- local
    integer(I4B) :: file_size
    !
    inquire (unit=locat, size=file_size)
    !
    if (expected_size /= file_size) then
      write (errmsg, '(a,i0,a,i0,a)') &
        'Unexpected file size for binary input array '// &
        trim(arrname)//'. Expected=', expected_size, &
        '/Found=', file_size, ' bytes.'
      call store_error(errmsg)
      call store_error_unit(locat)
    end if
  end subroutine check_binary_filesize

  !> @ brief Check the binary data size
  !!
  !!  Check the size of the binary data that will be read
  !!  relative to the unfilled elements in the array .
  !!
  !<
  function check_binary_size(nval, nvalt, arrsize, aname, locat) result(isok)
    ! -- dummy
    integer(I4B), intent(in) :: nval !< number of array
    integer(I4B), intent(in) :: nvalt !< current data index
    integer(I4B), intent(in) :: arrsize !< size of the array
    character(len=*), intent(in) :: aname !< name of array
    integer(I4B), intent(in) :: locat !< binary file unit
    !
    ! -- local variables
    logical(LGP) :: isok
    !
    ! -- initialize isok
    isok = .TRUE.
    !
    if (nvalt + nval > arrsize) then
      write (errmsg, '(a,i0,a,1x,a,1x,a,i0,a,1x,i0,3(1x,a))') &
        'The size of the data array calculated from the binary header (', &
        nval, ') will exceed the remainder of the', trim(adjustl(aname)), &
        'data array (', arrsize, ') array by', nvalt + nval - arrsize, &
        'elements. This is usually caused by incorrect assignment of', &
        '(m1,m2,m3) in the binary header. See the mf6io.pdf document', &
        'for information on assigning (m1,m2,m3).'
      call store_error(errmsg)
      call store_error_unit(locat)
      isok = .FALSE.
    end if
  end function check_binary_size

end module ArrayReadersModule
