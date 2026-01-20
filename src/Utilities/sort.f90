module SortModule
  use KindModule
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error
  use ConstantsModule, only: LINELENGTH
  implicit none

  private
  public :: qsort, selectn, unique_values

  interface qsort
    module procedure qsort_int1d, qsort_dbl1d
  end interface

  interface unique_values
    module procedure unique_values_int1d, unique_values_dbl1d
  end interface

contains

  !< @brief quick sort that also includes an index number
  !<
  subroutine qsort_int1d(indx, v, reverse)
    ! -- dummy arguments
    integer(I4B), dimension(:), intent(inout) :: indx
    integer(I4B), dimension(:), intent(inout) :: v
    logical, intent(in), optional :: reverse
    ! -- local variables
    logical :: lrev
    integer(I4B), parameter :: nn = 15
    integer(I4B), parameter :: nstack = 50
    integer(I4B) :: nsize
    integer(I4B) :: k
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: jstack
    integer(I4B) :: ileft
    integer(I4B) :: iright
    integer(I4B), dimension(nstack) :: istack
    integer(I4B) :: iidx
    integer(I4B) :: ia
    integer(I4B) :: a
    ! -- functions
    ! -- code
    !
    ! -- process optional dummy variables
    if (present(reverse)) then
      lrev = reverse
    else
      lrev = .FALSE.
    end if
    !
    ! -- initialize variables
    nsize = size(v)
    jstack = 0
    ileft = 1
    iright = nsize
    !
    ! -- perform quicksort
    do
      if (iright - ileft < nn) then
        do j = (ileft + 1), iright
          a = v(j)
          iidx = indx(j)
          do i = (j - 1), ileft, -1
            if (v(i) <= a) exit
            v(i + 1) = v(i)
            indx(i + 1) = indx(i)
          end do
          v(i + 1) = a
          indx(i + 1) = iidx
        end do
        if (jstack == 0) return
        iright = istack(jstack)
        ileft = istack(jstack - 1)
        jstack = jstack - 2
      else
        k = (ileft + iright) / 2
        call iswap(v(k), v(ileft + 1))
        call iswap(indx(k), indx(ileft + 1))
        if (v(ileft) > v(iright)) then
          call iswap(v(ileft), v(iright))
          call iswap(indx(ileft), indx(iright))
        end if
        if (v(ileft + 1) > v(iright)) then
          call iswap(v(ileft + 1), v(iright))
          call iswap(indx(ileft + 1), indx(iright))
        end if
        if (v(ileft) > v(ileft + 1)) then
          call iswap(v(ileft), v(ileft + 1))
          call iswap(indx(ileft), indx(ileft + 1))
        end if
        i = ileft + 1
        j = iright
        a = v(ileft + 1)
        ia = indx(ileft + 1)
        do
          do
            i = i + 1
            if (v(i) >= a) then
              exit
            end if
          end do
          do
            j = j - 1
            if (v(j) <= a) then
              exit
            end if
          end do
          if (j < i) then
            exit
          end if
          call iswap(v(i), v(j))
          call iswap(indx(i), indx(j))
        end do
        v(ileft + 1) = v(j)
        indx(ileft + 1) = indx(j)
        v(j) = a
        indx(j) = ia
        jstack = jstack + 2
        if (jstack > nstack) then
          write (errmsg, '(a,3(1x,a))') &
            'JSTACK > NSTACK IN SortModule::qsort'
          call store_error(errmsg, terminate=.TRUE.)
        end if
        if ((iright - i + 1) >= (j - 1)) then
          istack(jstack) = iright
          istack(jstack - 1) = i
          iright = j - 1
        else
          istack(jstack) = j - 1
          istack(jstack - 1) = ileft
          ileft = i
        end if
      end if
    end do
    !
    ! -- reverse order of the heap index
    if (lrev) then
      j = nsize
      do i = 1, nsize / 2
        call iswap(v(i), v(j))
        call iswap(indx(i), indx(j))
        j = j - 1
      end do
    end if
  end subroutine qsort_int1d

  !< @brief quick sort that also includes an index number
  !<
  subroutine qsort_dbl1d(indx, v, reverse)
    ! -- dummy arguments
    integer(I4B), dimension(:), intent(inout) :: indx
    real(DP), dimension(:), intent(inout) :: v
    logical, intent(in), optional :: reverse
    ! -- local variables
    logical :: lrev
    integer(I4B), parameter :: nn = 15
    integer(I4B), parameter :: nstack = 50
    integer(I4B) :: nsize
    integer(I4B) :: k
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: jstack
    integer(I4B) :: ileft
    integer(I4B) :: iright
    integer(I4B), dimension(nstack) :: istack
    integer(I4B) :: iidx
    integer(I4B) :: ia
    real(DP) :: a
    ! -- functions
    ! -- code
    !
    ! -- process optional dummy variables
    if (present(reverse)) then
      lrev = reverse
    else
      lrev = .FALSE.
    end if
    !
    ! -- initialize variables
    nsize = size(v)
    jstack = 0
    ileft = 1
    iright = nsize
    !
    ! -- perform quicksort
    do
      if (iright - ileft < nn) then
        do j = (ileft + 1), iright
          a = v(j)
          iidx = indx(j)
          do i = (j - 1), ileft, -1
            if (v(i) <= a) exit
            v(i + 1) = v(i)
            indx(i + 1) = indx(i)
          end do
          v(i + 1) = a
          indx(i + 1) = iidx
        end do
        if (jstack == 0) return
        iright = istack(jstack)
        ileft = istack(jstack - 1)
        jstack = jstack - 2
      else
        k = (ileft + iright) / 2
        call rswap(v(k), v(ileft + 1))
        call iswap(indx(k), indx(ileft + 1))
        if (v(ileft) > v(iright)) then
          call rswap(v(ileft), v(iright))
          call iswap(indx(ileft), indx(iright))
        end if
        if (v(ileft + 1) > v(iright)) then
          call rswap(v(ileft + 1), v(iright))
          call iswap(indx(ileft + 1), indx(iright))
        end if
        if (v(ileft) > v(ileft + 1)) then
          call rswap(v(ileft), v(ileft + 1))
          call iswap(indx(ileft), indx(ileft + 1))
        end if
        i = ileft + 1
        j = iright
        a = v(ileft + 1)
        ia = indx(ileft + 1)
        do
          do
            i = i + 1
            if (v(i) >= a) then
              exit
            end if
          end do
          do
            j = j - 1
            if (v(j) <= a) then
              exit
            end if
          end do
          if (j < i) then
            exit
          end if
          call rswap(v(i), v(j))
          call iswap(indx(i), indx(j))
        end do
        v(ileft + 1) = v(j)
        indx(ileft + 1) = indx(j)
        v(j) = a
        indx(j) = ia
        jstack = jstack + 2
        if (jstack > nstack) then
          write (errmsg, '(a,3(1x,a))') &
            'JSTACK > NSTACK IN SortModule::qsort'
          call store_error(errmsg, terminate=.TRUE.)
        end if
        if ((iright - i + 1) >= (j - 1)) then
          istack(jstack) = iright
          istack(jstack - 1) = i
          iright = j - 1
        else
          istack(jstack) = j - 1
          istack(jstack - 1) = ileft
          ileft = i
        end if
      end if
    end do
    !
    ! -- reverse order of the heap index
    if (lrev) then
      j = nsize
      do i = 1, nsize / 2
        call rswap(v(i), v(j))
        call iswap(indx(i), indx(j))
        j = j - 1
      end do
    end if
  end subroutine qsort_dbl1d

  subroutine unique_values_int1d(a, b)
    ! - dummy arguments
    integer(I4B), dimension(:), allocatable, intent(in) :: a
    integer(I4B), dimension(:), allocatable, intent(inout) :: b
    ! -- local variables
    integer(I4B) :: count
    integer(I4B) :: n
    integer(I4B), dimension(:), allocatable :: indxarr
    integer(I4B), dimension(:), allocatable :: tarr
    ! -- functions
    ! -- code
    !
    ! -- allocate tarr and create idxarr
    allocate (tarr(size(a)))
    allocate (indxarr(size(a)))
    !
    ! -- fill tarr with a and create index
    do n = 1, size(a)
      tarr(n) = a(n)
      indxarr(n) = n
    end do
    !
    ! -- sort a in increasing order
    call qsort(indxarr, tarr, reverse=.TRUE.)
    !
    ! -- determine the number of unique values
    count = 1
    do n = 2, size(tarr)
      if (tarr(n) > tarr(n - 1)) count = count + 1
    end do
    !
    ! -- allocate b for unique values
    if (allocated(b)) then
      deallocate (b)
    end if
    allocate (b(count))
    !
    ! -- fill b with unique values
    b(1) = tarr(1)
    count = 1
    do n = 2, size(a)
      if (tarr(n) > b(count)) then
        count = count + 1
        b(count) = tarr(n)
      end if
    end do
    !
    ! -- allocate tarr and create idxarr
    deallocate (tarr)
    deallocate (indxarr)
  end subroutine unique_values_int1d

  subroutine unique_values_dbl1d(a, b)
    ! - dummy arguments
    real(DP), dimension(:), allocatable, intent(in) :: a
    real(DP), dimension(:), allocatable, intent(inout) :: b
    ! -- local variables
    integer(I4B) :: count
    integer(I4B) :: n
    integer(I4B), dimension(:), allocatable :: indxarr
    real(DP), dimension(:), allocatable :: tarr
    ! -- functions
    ! -- code
    !
    ! -- allocate tarr and create idxarr
    allocate (tarr(size(a)))
    allocate (indxarr(size(a)))
    !
    ! -- fill tarr with a and create index
    do n = 1, size(a)
      tarr(n) = a(n)
      indxarr(n) = n
    end do
    !
    ! -- sort a in increasing order
    call qsort(indxarr, tarr, reverse=.TRUE.)
    !
    ! -- determine the number of unique values
    count = 1
    do n = 2, size(tarr)
      if (tarr(n) > tarr(n - 1)) count = count + 1
    end do
    !
    ! -- allocate b for unique values
    if (allocated(b)) then
      deallocate (b)
    end if
    allocate (b(count))
    !
    ! -- fill b with unique values
    b(1) = tarr(1)
    count = 1
    do n = 2, size(a)
      if (tarr(n) > b(count)) then
        count = count + 1
        b(count) = tarr(n)
      end if
    end do
    !
    ! -- allocate tarr and create idxarr
    deallocate (tarr)
    deallocate (indxarr)
  end subroutine unique_values_dbl1d

  !< @brief heap selection
  !<
  subroutine selectn(indx, v, reverse)
    ! -- dummy arguments
    integer(I4B), dimension(:), intent(inout) :: indx
    real(DP), dimension(:), intent(inout) :: v
    logical, intent(in), optional :: reverse
    ! -- local variables
    logical :: lrev
    integer(I4B) :: nsizei
    integer(I4B) :: nsizev
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: k
    integer(I4B) :: n
    !integer(I4B) :: iidx
    real(DP), dimension(:), allocatable :: vv
    ! -- functions
    ! -- code
    !
    ! -- process optional dummy variables
    if (present(reverse)) then
      lrev = reverse
    else
      lrev = .FALSE.
    end if
    !
    ! -- initialize heap
    nsizev = size(v)
    nsizei = min(nsizev, size(indx))
    allocate (vv(nsizei))
    !
    ! -- initialize heap index (indx) and heap (vv)
    do n = 1, nsizei
      vv(n) = v(n)
      indx(n) = n
    end do
    !
    ! -- initial sort
    call qsort(indx, vv)
    !
    ! -- evaluate the remaining elements in v
    do i = nsizei + 1, nsizev
      !
      ! -- put the current value on the heap
      if (v(i) > vv(1)) then
        vv(1) = v(i)
        indx(1) = i
        j = 1
        do
          k = 2 * j
          if (k > nsizei) then
            exit
          end if
          if (k /= nsizei) then
            if (vv(k) > vv(k + 1)) then
              k = k + 1
            end if
          end if
          if (vv(j) <= vv(k)) then
            exit
          end if
          call rswap(vv(k), vv(j))
          call iswap(indx(k), indx(j))
          j = k
        end do
      end if
    end do
    !
    ! -- final sort
    call qsort(indx, vv)
    !
    ! -- reverse order of the heap index
    if (lrev) then
      j = nsizei
      do i = 1, nsizei / 2
        call iswap(indx(i), indx(j))
        j = j - 1
      end do
    end if
  end subroutine selectn

  subroutine rswap(a, b)
    ! -- dummy arguments
    real(DP), intent(inout) :: a
    real(DP), intent(inout) :: b
    ! -- local variables
    real(DP) :: d
    ! -- functions
    ! -- code
    d = a
    a = b
    b = d
  end subroutine rswap

  subroutine iswap(ia, ib)
    ! -- dummy arguments
    integer(I4B), intent(inout) :: ia
    integer(I4B), intent(inout) :: ib
    ! -- local variables
    integer(I4B) :: id
    ! -- functions
    ! -- code
    id = ia
    ia = ib
    ib = id
  end subroutine iswap

end module SortModule
