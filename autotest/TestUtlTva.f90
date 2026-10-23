module TestUtlTva
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use KindModule, only: I4B, DP
  use ConstantsModule, only: LENAUXNAME, DZERO
  use CharacterStringModule, only: CharacterStringType
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryHelperModule, only: create_mem_path
  use BaseDisModule, only: DisBaseType
  use UtlTvaModule, only: UtlTvaType
  implicit none
  private
  public :: collect_utltva

  ! a three-layer, two-column grid: the point is that TVA must address every
  ! node, not just the first layer
  integer(I4B), parameter :: NLAY = 3
  integer(I4B), parameter :: NCPL = 2
  integer(I4B), parameter :: NODES = NLAY * NCPL
  integer(I4B), parameter :: NVAR = 2

contains

  subroutine collect_utltva(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
                new_unittest("has", test_has), &
                new_unittest("fill_full_grid", test_fill_full_grid), &
                new_unittest("fill_before_period", test_fill_before_period) &
                ]
  end subroutine collect_utltva

  !> @brief Minimal discretization: TVA only uses nodes and get_nodeuser, and
  !! get_nodeuser is the identity when nodes is not less than nodesuser.
  !<
  subroutine make_dis(dis)
    class(DisBaseType), pointer, intent(out) :: dis
    allocate (dis)
    allocate (dis%nodes)
    allocate (dis%nodesuser)
    dis%nodes = NODES
    dis%nodesuser = NODES
  end subroutine make_dis

  subroutine free_dis(dis)
    class(DisBaseType), pointer, intent(inout) :: dis
    deallocate (dis%nodes)
    deallocate (dis%nodesuser)
    deallocate (dis)
  end subroutine free_dis

  !> @brief Stand in for the Input Data Model by putting AUXILIARY, and
  !! optionally AUXVAR, where the loader would have put them.
  !<
  subroutine make_input(mempath, with_auxvar)
    character(len=*), intent(in) :: mempath
    logical, intent(in) :: with_auxvar
    type(CharacterStringType), dimension(:), pointer, contiguous :: aux
    real(DP), dimension(:, :), pointer, contiguous :: auxvar
    integer(I4B) :: iv, n
    call mem_allocate(aux, LENAUXNAME, NVAR, 'AUXILIARY', mempath)
    aux(1) = 'SALTWATER_HEAD'
    aux(2) = 'OTHER'
    if (with_auxvar) then
      call mem_allocate(auxvar, NVAR, NODES, 'AUXVAR', mempath)
      do iv = 1, NVAR
        do n = 1, NODES
          auxvar(iv, n) = real(100 * iv + n, DP)
        end do
      end do
    end if
  end subroutine make_input

  subroutine free_input(mempath, with_auxvar)
    character(len=*), intent(in) :: mempath
    logical, intent(in) :: with_auxvar
    type(CharacterStringType), dimension(:), pointer, contiguous :: aux
    real(DP), dimension(:, :), pointer, contiguous :: auxvar
    call mem_deallocate(aux, 'AUXILIARY', mempath)
    if (with_auxvar) call mem_deallocate(auxvar, 'AUXVAR', mempath)
  end subroutine free_input

  subroutine test_has(error)
    type(error_type), allocatable, intent(out) :: error
    class(DisBaseType), pointer :: dis
    type(UtlTvaType) :: tva
    character(len=*), parameter :: IN = 'TSTA/TVAIN'

    call make_dis(dis)
    call make_input(IN, .true.)
    call tva%initialize(dis, NODES, IN, 0, 'TSTA', 'PKGA')

    call check(error, tva%has('SALTWATER_HEAD'), &
               "a declared auxiliary name should be found")
    if (allocated(error)) return
    call check(error,.not. tva%has('NOT_THERE'), &
               "an undeclared name should not be found")

    call tva%tva_da()
    call free_input(IN, .true.)
    call free_dis(dis)
  end subroutine test_has

  !> @brief fill must copy values for every node in the grid, not only the
  !! first layer. This is the regression that motivated full-grid support.
  !<
  subroutine test_fill_full_grid(error)
    type(error_type), allocatable, intent(out) :: error
    class(DisBaseType), pointer :: dis
    type(UtlTvaType) :: tva
    real(DP), dimension(NODES) :: vec
    integer(I4B) :: n
    character(len=*), parameter :: IN = 'TSTB/TVAIN'

    call make_dis(dis)
    call make_input(IN, .true.)
    call tva%initialize(dis, NODES, IN, 0, 'TSTB', 'PKGB')

    vec = DZERO
    call tva%fill(vec, 'SALTWATER_HEAD')
    do n = 1, NODES
      call check(error, vec(n), real(100 + n, DP), &
                 "fill should copy the first variable for every node")
      if (allocated(error)) return
    end do

    vec = DZERO
    call tva%fill(vec, 'OTHER')
    do n = 1, NODES
      call check(error, vec(n), real(200 + n, DP), &
                 "fill should select by auxiliary name")
      if (allocated(error)) return
    end do

    call tva%tva_da()
    call free_input(IN, .true.)
    call free_dis(dis)
  end subroutine test_fill_full_grid

  !> @brief Before any period block supplies values the loaded array is
  !! absent, and fill must leave the caller's vector untouched.
  !<
  subroutine test_fill_before_period(error)
    type(error_type), allocatable, intent(out) :: error
    class(DisBaseType), pointer :: dis
    type(UtlTvaType) :: tva
    real(DP), dimension(NODES) :: vec
    integer(I4B) :: n
    character(len=*), parameter :: IN = 'TSTC/TVAIN'

    call make_dis(dis)
    call make_input(IN, .false.)
    call tva%initialize(dis, NODES, IN, 0, 'TSTC', 'PKGC')

    vec = -1.0_DP
    call tva%fill(vec, 'SALTWATER_HEAD')
    do n = 1, NODES
      call check(error, vec(n), -1.0_DP, &
                 "fill should not modify the vector before a period block")
      if (allocated(error)) return
    end do

    call tva%tva_da()
    call free_input(IN, .false.)
    call free_dis(dis)
  end subroutine test_fill_before_period

end module TestUtlTva
