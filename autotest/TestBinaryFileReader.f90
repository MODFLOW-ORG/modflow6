module TestBinaryFileReader
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use KindModule, only: I4B, I8B, DP, LGP
  use BudgetFileReaderModule, only: BudgetFileReaderType, BudgetFileHeaderType
  use HeadFileReaderModule, only: HeadFileReaderType, HeadFileHeaderType
  implicit none
  private
  public :: collect_binaryfilereader

  !> First record starts just below the 2 GiB signed 32-bit boundary so that
  !! it straddles the boundary and the second record starts beyond it. Bytes
  !! before this position are never written, so on most file systems the test
  !! file is sparse and uses almost no disk space.
  integer(I8B), parameter :: POS0 = 2_I8B**31 - 100_I8B

contains

  subroutine collect_binaryfilereader(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("budget_file_beyond_2gib", &
                             test_budget_file_beyond_2gib), &
                new_unittest("head_file_beyond_2gib", &
                             test_head_file_beyond_2gib) &
                ]
  end subroutine collect_binaryfilereader

  !> @brief Read budget records located across and beyond the 2 GiB offset
  subroutine test_budget_file_beyond_2gib(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: fname = 'test_bfr_2gib.cbc'
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
    open (newunit=iu, file=fname, access='stream', form='unformatted', &
          status='replace', action='write')
    write (iu, pos=POS0) 1, 1, '    FLOW-JA-FACE', nja, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    write (iu) flowja
    inquire (unit=iu, pos=pos1)
    write (iu) 1, 1, '          STO-SS', nflow, 1, -1
    write (iu) 1, 1.0_DP, 1.0_DP, 1.0_DP
    write (iu) flow
    close (iu)

    ! position the reader at the first record
    open (newunit=iu, file=fname, access='stream', form='unformatted', &
          status='old', action='read')
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

100 close (iu, status='delete')
  end subroutine test_budget_file_beyond_2gib

  !> @brief Read head records located across and beyond the 2 GiB offset
  subroutine test_head_file_beyond_2gib(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: fname = 'test_hfr_2gib.hds'
    integer(I4B), parameter :: ncol = 10, nrow = 5
    type(HeadFileReaderType) :: hfr
    real(DP) :: head1(ncol * nrow), head2(ncol * nrow)
    integer(I8B) :: pos1
    integer(I4B) :: iu, i
    logical(LGP) :: success

    head1 = [(real(i, DP), i=1, ncol * nrow)]
    head2 = -head1

    ! write two head records starting at POS0
    open (newunit=iu, file=fname, access='stream', form='unformatted', &
          status='replace', action='write')
    write (iu, pos=POS0) 1, 1, 1.0_DP, 1.0_DP, '            HEAD', ncol, nrow, 1
    write (iu) head1
    inquire (unit=iu, pos=pos1)
    write (iu) 2, 1, 2.0_DP, 2.0_DP, '            HEAD', ncol, nrow, 1
    write (iu) head2
    close (iu)

    ! position the reader at the first record
    open (newunit=iu, file=fname, access='stream', form='unformatted', &
          status='old', action='read')
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

100 close (iu, status='delete')
  end subroutine test_head_file_beyond_2gib

end module TestBinaryFileReader
