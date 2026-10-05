module TestGridFileReader
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use KindModule, only: I4B, I8B, DP
  use GridFileReaderModule, only: GridFileReaderType
  implicit none
  private
  public :: collect_gridfilereader

  integer(I4B), parameter :: LENHDR = 50
  integer(I4B), parameter :: LENTXT = 100

contains

  subroutine collect_gridfilereader(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("grid_file_beyond_2gib", &
                             test_grid_file_beyond_2gib) &
                ]
  end subroutine collect_gridfilereader

  !> @brief Write a header or variable definition line to a grid file
  subroutine write_line(iu, str, n)
    integer(I4B), intent(in) :: iu
    character(len=*), intent(in) :: str
    integer(I4B), intent(in) :: n
    character(len=n) :: line

    line = str
    line(n:n) = new_line('a')
    write (iu) line
  end subroutine write_line

  !> @brief Read variables located beyond the 2 GiB offset
  !!
  !! The reader indexes whatever the header declares, so declare a fake
  !! 2 GiB PAD array to push the other variables past 2 GiB. PAD is never
  !! written, so the file is sparse where the file system supports it.
  !<
  subroutine test_grid_file_beyond_2gib(error)
    type(error_type), allocatable, intent(out) :: error
    integer(I4B), parameter :: nodes = 5
    integer(I4B), parameter :: npad = 2**29
    character(len=*), parameter :: crs = 'EPSG:26916'
    type(GridFileReaderType) :: gfr
    integer(I4B) :: idomain(nodes), idomain_read(nodes)
    real(DP) :: botm(nodes)
    real(DP), allocatable :: botm_read(:)
    character(len=:), allocatable :: crs_read
    integer(I8B) :: pos
    integer(I4B) :: iu, i
    character(len=LENTXT) :: txt

    idomain = [(i, i=1, nodes)]
    botm = [(-real(i, DP), i=1, nodes)]

    ! header
    open (newunit=iu, access='stream', form='unformatted', status='scratch')
    call write_line(iu, 'GRID DISU', LENHDR)
    call write_line(iu, 'VERSION 2', LENHDR)
    call write_line(iu, 'NTXT 7', LENHDR)
    write (txt, '(a, i0)') 'LENTXT ', LENTXT
    call write_line(iu, txt, LENHDR)
    write (txt, '(a, i0)') 'NODES INTEGER NDIM 0 # ', nodes
    call write_line(iu, txt, LENTXT)
    write (txt, '(a, i0)') 'PAD INTEGER NDIM 1 ', npad
    call write_line(iu, txt, LENTXT)
    write (txt, '(a, i0)') 'IDOMAIN INTEGER NDIM 1 ', nodes
    call write_line(iu, txt, LENTXT)
    write (txt, '(a, i0)') 'BOTM DOUBLE NDIM 1 ', nodes
    call write_line(iu, txt, LENTXT)
    call write_line(iu, 'ANGROT DOUBLE NDIM 0 # 30.0', LENTXT)
    write (txt, '(a, i0)') 'CRS CHARACTER NDIM 1 ', len(crs)
    call write_line(iu, txt, LENTXT)
    call write_line(iu, 'NCPL INTEGER NDIM 0 # 7', LENTXT)

    ! data, skipping over the padding
    write (iu) nodes
    inquire (unit=iu, pos=pos)
    pos = pos + int(npad, I8B) * 4
    write (iu, pos=pos) idomain
    write (iu) botm
    write (iu) 30.0_DP
    write (iu) crs
    write (iu) 7
    rewind (iu)

    call gfr%initialize(iu)

    call check(error, gfr%read_int('NODES') == nodes, 'wrong NODES')
    if (allocated(error)) goto 100
    call gfr%read_int_1d_into('IDOMAIN', idomain_read)
    call check(error, all(idomain_read == idomain), 'wrong IDOMAIN')
    if (allocated(error)) goto 100
    botm_read = gfr%read_dbl_1d('BOTM')
    call check(error, all(botm_read == botm), 'wrong BOTM')
    if (allocated(error)) goto 100
    call check(error, gfr%read_dbl('ANGROT') == 30.0_DP, 'wrong ANGROT')
    if (allocated(error)) goto 100
    crs_read = gfr%read_charstr('CRS')
    call check(error, crs_read == crs, 'wrong CRS')
    if (allocated(error)) goto 100
    call check(error, gfr%read_int('NCPL') == 7, 'wrong NCPL')

100 call gfr%finalize()
  end subroutine test_grid_file_beyond_2gib

end module TestGridFileReader
