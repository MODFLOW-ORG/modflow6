module TestBinaryFileReader
  use testdrive, only: error_type, unittest_type, new_unittest, check, &
                       skip_test
  use KindModule, only: I4B, I8B, DP, LGP
  use BudgetFileReaderModule, only: BudgetFileReaderType, BudgetFileHeaderType
  use HeadFileReaderModule, only: HeadFileReaderType, HeadFileHeaderType
  implicit none
  private
  public :: collect_binaryfilereader

  !> First record starts just below the 2 GiB signed 32-bit boundary so that
  !! it straddles the boundary and the second record starts beyond it. Bytes
  !! before this position are never written, so on most file systems the test
  !! file is sparse and uses almost no disk space. Scratch files give each run
  !! a unique file, since meson may run the suite in parallel processes.
  integer(I8B), parameter :: POS0 = 2_I8B**31 - 100_I8B

  !> Set this environment variable to run tests that read 2 GiB or more
  character(len=*), parameter :: LARGE_TESTS_VAR = 'MF6_LARGE_TESTS'

contains

  subroutine collect_binaryfilereader(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("budget_file_beyond_2gib", &
                             test_budget_file_beyond_2gib), &
                new_unittest("head_file_beyond_2gib", &
                             test_head_file_beyond_2gib), &
                new_unittest("budget_file_index", &
                             test_budget_file_index), &
                new_unittest("budget_file_index_beyond_2gib", &
                             test_budget_file_index_beyond_2gib) &
                ]
  end subroutine collect_binaryfilereader

  !> @brief Read budget records located across and beyond the 2 GiB offset
  subroutine test_budget_file_beyond_2gib(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: nja = 50
    integer(I4B), parameter :: nflow = 3
    type(BudgetFileReaderType) :: bfr
    real(DP) :: flowja(nja), flow(nflow)
    integer(I8B) :: pos1
    integer(I4B) :: iu, i
    logical(LGP) :: success

    flowja = [(real(i, DP), i=1, nja)]
    flow = [(-real(i, DP), i=1, nflow)]

    ! write two imeth=1 records starting at POS0
    open (newunit=iu, access='stream', form='unformatted', status='scratch')
    write (iu, pos=POS0) 1, 1, '    FLOW-JA-FACE', nja, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    write (iu) flowja
    inquire (unit=iu, pos=pos1)
    write (iu) 1, 1, '          STO-SS', nflow, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    write (iu) flow

    ! position the reader at the first record
    bfr%inunit = iu
    call bfr%rewind()
    read (iu, pos=POS0)

    ! first record straddles 2 GiB
    call bfr%read_record(success)
    call check(error, success, 'failed to read first record')
    if (allocated(error)) goto 100
    call check(error, bfr%header%pos == POS0, 'wrong first record position')
    if (allocated(error)) goto 100
    call check(error, all(bfr%flowja == flowja), 'wrong FLOW-JA-FACE values')
    if (allocated(error)) goto 100
    call check(error,.not. bfr%endoffile, 'second record not found')
    if (allocated(error)) goto 100

    ! second record lies entirely beyond 2 GiB
    call bfr%read_record(success)
    call check(error, success, 'failed to read second record')
    if (allocated(error)) goto 100
    call check(error, bfr%header%pos == pos1, 'wrong second record position')
    if (allocated(error)) goto 100
    select type (h => bfr%header)
    type is (BudgetFileHeaderType)
      call check(error, h%budtxt == '          STO-SS', 'wrong budget text')
      if (allocated(error)) goto 100
      call check(error, h%imeth == 1, 'wrong method code')
      if (allocated(error)) goto 100
    end select
    call check(error, all(bfr%flow == flow), 'wrong STO-SS values')
    if (allocated(error)) goto 100
    call check(error, bfr%endoffile, 'end of file not detected')

100 close (iu)
  end subroutine test_budget_file_beyond_2gib

  !> @brief Read head records located across and beyond the 2 GiB offset
  subroutine test_head_file_beyond_2gib(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: ncol = 10, nrow = 5
    type(HeadFileReaderType) :: hfr
    real(DP) :: head1(ncol * nrow), head2(ncol * nrow)
    integer(I8B) :: pos1
    integer(I4B) :: iu, i
    logical(LGP) :: success

    head1 = [(real(i, DP), i=1, ncol * nrow)]
    head2 = -head1

    ! write two head records starting at POS0
    open (newunit=iu, access='stream', form='unformatted', status='scratch')
    write (iu, pos=POS0) 1, 1, 1.0_DP, 1.0_DP, '            HEAD', ncol, nrow, 1
    write (iu) head1
    inquire (unit=iu, pos=pos1)
    write (iu) 2, 1, 2.0_DP, 2.0_DP, '            HEAD', ncol, nrow, 1
    write (iu) head2

    ! position the reader at the first record
    hfr%inunit = iu
    call hfr%rewind()
    read (iu, pos=POS0)

    ! first record straddles 2 GiB
    call hfr%read_record(success)
    call check(error, success, 'failed to read first record')
    if (allocated(error)) goto 100
    call check(error, hfr%header%pos == POS0, 'wrong first record position')
    if (allocated(error)) goto 100
    call check(error, all(hfr%head == head1), 'wrong first head values')
    if (allocated(error)) goto 100
    call check(error,.not. hfr%endoffile, 'second record not found')
    if (allocated(error)) goto 100

    ! second record lies entirely beyond 2 GiB
    call hfr%read_record(success)
    call check(error, success, 'failed to read second record')
    if (allocated(error)) goto 100
    call check(error, hfr%header%pos == pos1, 'wrong second record position')
    if (allocated(error)) goto 100
    call check(error, hfr%header%kstp == 2, 'wrong second record kstp')
    if (allocated(error)) goto 100
    call check(error, all(hfr%head == head2), 'wrong second head values')
    if (allocated(error)) goto 100
    call check(error, hfr%endoffile, 'end of file not detected')

100 close (iu)
  end subroutine test_head_file_beyond_2gib

  !> @brief Index budget records and seek to them by index
  subroutine test_budget_file_index(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: nrec = 3, nval = 4
    character(len=16), parameter :: budtxt(nrec) = [ &
                                    '          STO-SS', &
                                    '          STO-SY', &
                                    '             WEL']
    type(BudgetFileReaderType) :: bfr
    real(DP) :: flow(nval, nrec)
    integer(I8B) :: expected(nrec)
    integer(I4B) :: iu, i, k
    logical(LGP) :: success

    flow = reshape([(real(i, DP), i=1, nval * nrec)], [nval, nrec])

    ! write records from the start of the file, saving their positions
    open (newunit=iu, access='stream', form='unformatted', status='scratch')
    do k = 1, nrec
      inquire (unit=iu, pos=expected(k))
      write (iu) 1, 1, budtxt(k), nval, 1, -1
      write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
      write (iu) flow(:, k)
    end do

    bfr%inunit = iu
    call bfr%build_index()
    call check(error, bfr%indexed, 'file not indexed')
    if (allocated(error)) goto 100
    call check(error, bfr%nrecords == nrec, 'wrong number of records')
    if (allocated(error)) goto 100
    call check(error, all(bfr%record_positions == expected), &
               'wrong record positions')
    if (allocated(error)) goto 100

    ! seek to the second record and read it
    call bfr%seek_to_index(2)
    call check(error,.not. bfr%endoffile, 'end of file after seek')
    if (allocated(error)) goto 100
    call bfr%read_record(success)
    call check(error, success, 'failed to read second record')
    if (allocated(error)) goto 100
    call check(error, bfr%header%pos == expected(2), &
               'wrong second record position')
    if (allocated(error)) goto 100
    select type (h => bfr%header)
    type is (BudgetFileHeaderType)
      call check(error, h%budtxt == budtxt(2), 'wrong budget text')
      if (allocated(error)) goto 100
    end select
    call check(error, all(bfr%flow == flow(:, 2)), 'wrong second record values')
    if (allocated(error)) goto 100

    ! seeking past the last record signals end of file
    call bfr%seek_to_index(nrec + 1)
    call check(error, bfr%endoffile, 'end of file not signaled')

100 close (iu)
  end subroutine test_budget_file_index

  !> @brief Index budget records located beyond the 2 GiB offset
  !!
  !! The first record holds 2 GiB of FLOW-JA-FACE data. Only its header is
  !! written, so the data region is sparse on most file systems, but indexing
  !! reads it twice and allocates it in memory. Skipped unless the
  !! MF6_LARGE_TESTS environment variable is set.
  subroutine test_budget_file_index_beyond_2gib(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: nja = 2**28
    integer(I4B), parameter :: nrec = 3, nval = 4
    character(len=16), parameter :: budtxt(nrec) = [ &
                                    '    FLOW-JA-FACE', &
                                    '          STO-SS', &
                                    '             WEL']
    type(BudgetFileReaderType) :: bfr
    real(DP) :: flow(nval, 2:nrec)
    integer(I8B) :: expected(nrec), pos
    integer(I4B) :: iu, i, k, envlen
    logical(LGP) :: success

    call get_environment_variable(LARGE_TESTS_VAR, length=envlen)
    if (envlen == 0) then
      call skip_test(error, 'set '//LARGE_TESTS_VAR//' to run')
      return
    end if

    flow = reshape([(real(i, DP), i=1, nval * (nrec - 1))], [nval, nrec - 1])

    ! write the first record's header, then skip past its 2 GiB of data
    open (newunit=iu, access='stream', form='unformatted', status='scratch')
    inquire (unit=iu, pos=expected(1))
    write (iu) 1, 1, budtxt(1), nja, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    inquire (unit=iu, pos=pos)
    pos = pos + int(nja, I8B) * int(storage_size(1.0_DP) / 8, I8B)

    ! write the remaining records after the skipped data
    expected(2) = pos
    write (iu, pos=pos) 1, 1, budtxt(2), nval, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    write (iu) flow(:, 2)
    inquire (unit=iu, pos=expected(3))
    write (iu) 1, 1, budtxt(3), nval, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    write (iu) flow(:, 3)

    bfr%inunit = iu
    call bfr%build_index()
    call check(error, bfr%nrecords == nrec, 'wrong number of records')
    if (allocated(error)) goto 100
    call check(error, all(bfr%record_positions == expected), &
               'wrong record positions')
    if (allocated(error)) goto 100
    call check(error, all(bfr%record_positions(2:) > 2_I8B**31), &
               'records not beyond 2 GiB')
    if (allocated(error)) goto 100

    ! seek to the last record and read it
    call bfr%seek_to_index(nrec)
    call bfr%read_record(success)
    call check(error, success, 'failed to read last record')
    if (allocated(error)) goto 100
    call check(error, bfr%header%pos == expected(nrec), &
               'wrong last record position')
    if (allocated(error)) goto 100
    select type (h => bfr%header)
    type is (BudgetFileHeaderType)
      call check(error, h%budtxt == budtxt(nrec), 'wrong budget text')
      if (allocated(error)) goto 100
    end select
    call check(error, all(bfr%flow == flow(:, nrec)), 'wrong last record values')
    if (allocated(error)) goto 100
    call check(error, bfr%endoffile, 'end of file not detected')

100 close (iu)
  end subroutine test_budget_file_index_beyond_2gib

end module TestBinaryFileReader
