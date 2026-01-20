MODULE IMSReorderingModule
  use KindModule, only: DP, I4B
  private
  public :: ims_odrv
contains

  subroutine ims_odrv(n, nja, nsp, ia, ja, p, ip, isp, flag)
    !
    !                                                                3/12/82
    !***********************************************************************
    !  odrv -- driver for sparse matrix reordering routines
    !***********************************************************************
    !
    !  description
    !
    !    odrv finds a minimum degree ordering of the rows and columns
    !    of a matrix m stored in (ia,ja,a) format (see below).  for the
    !    reordered matrix, the work and storage required to perform
    !    gaussian elimination is (usually) significantly less.
    !
    !    note.. odrv and its subordinate routines have been modified to
    !    compute orderings for general matrices, not necessarily having any
    !    symmetry.  the minimum degree ordering is computed for the
    !    structure of the symmetric matrix  m + m-transpose.
    !    modifications to the original odrv module have been made in
    !    the coding in subroutine mdi, and in the initial comments in
    !    subroutines odrv and md.
    !
    !    if only the nonzero entries in the upper triangle of m are being
    !    stored, then odrv symmetrically reorders (ia,ja,a), (optionally)
    !    with the diagonal entries placed first in each row.  this is to
    !    ensure that if m(i,j) will be in the upper triangle of m with
    !    respect to the new ordering, then m(i,j) is stored in row i (and
    !    thus m(j,i) is not stored),  whereas if m(i,j) will be in the
    !    strict lower triangle of m, then m(j,i) is stored in row j (and
    !    thus m(i,j) is not stored).
    !
    !
    !  storage of sparse matrices
    !
    !    the nonzero entries of the matrix m are stored row-by-row in the
    !    array a.  to identify the individual nonzero entries in each row,
    !    we need to know in which column each entry lies.  these column
    !    indices are stored in the array ja.  i.e., if  a(k) = m(i,j),  then
    !    ja(k) = j.  to identify the individual rows, we need to know where
    !    each row starts.  these row pointers are stored in the array ia.
    !    i.e., if m(i,j) is the first nonzero entry (stored) in the i-th row
    !    and  a(k) = m(i,j),  then  ia(i) = k.  moreover, ia(n+1) points to
    !    the first location following the last element in the last row.
    !    thus, the number of entries in the i-th row is  ia(i+1) - ia(i),
    !    the nonzero entries in the i-th row are stored consecutively in
    !
    !            a(ia(i)),  a(ia(i)+1),  ..., a(ia(i+1)-1),
    !
    !    and the corresponding column indices are stored consecutively in
    !
    !            ja(ia(i)), ja(ia(i)+1), ..., ja(ia(i+1)-1).
    !
    !    since the coefficient matrix is symmetric, only the nonzero entries
    !    in the upper triangle need be stored.  for example, the matrix
    !
    !             ( 1  0  2  3  0 )
    !             ( 0  4  0  0  0 )
    !         m = ( 2  0  5  6  0 )
    !             ( 3  0  6  7  8 )
    !             ( 0  0  0  8  9 )
    !
    !    could be stored as
    !
    !            - 1  2  3  4  5  6  7  8  9 10 11 12 13
    !         ---+--------------------------------------
    !         ia - 1  4  5  8 12 14
    !         ja - 1  3  4  2  1  3  4  1  3  4  5  4  5
    !          a - 1  2  3  4  2  5  6  3  6  7  8  8  9
    !
    !    or (symmetrically) as
    !
    !            - 1  2  3  4  5  6  7  8  9
    !         ---+--------------------------
    !         ia - 1  4  5  7  9 10
    !         ja - 1  3  4  2  3  4  4  5  5
    !          a - 1  2  3  4  5  6  7  8  9          .
    !
    !
    !  parameters
    !
    !    n    - order of the matrix
    !
    !    nja  - number of nonzeroes in the matrix
    !
    !    nsp  - declared dimension of the one-dimensional array isp.  nsp
    !           must be at least  3n+4k,  where k is the number of nonzeroes
    !           in the strict upper triangle of m
    !
    !    ia   - integer one-dimensional array containing pointers to delimit
    !           rows in ja and a.  dimension = n+1
    !
    !    ja   - integer one-dimensional array containing the column indices
    !           corresponding to the elements of a.  dimension = number of
    !           nonzero entries in (the upper triangle of) m
    !
    !    a    - real one-dimensional array containing the nonzero entries in
    !           (the upper triangle of) m, stored by rows.  dimension =
    !           number of nonzero entries in (the upper triangle of) m
    !
    !    p    - integer one-dimensional array used to return the permutation
    !           of the rows and columns of m corresponding to the minimum
    !           degree ordering.  dimension = n
    !
    !    ip   - integer one-dimensional array used to return the inverse of
    !           the permutation returned in p.  dimension = n
    !
    !    isp  - integer one-dimensional array used for working storage.
    !           dimension = nsp
    !
    !    path - integer path specification.  values and their meanings are -
    !             1  find minimum degree ordering only
    !             2  find minimum degree ordering and reorder symmetrically
    !                  stored matrix (used when only the nonzero entries in
    !                  the upper triangle of m are being stored)
    !             3  reorder symmetrically stored matrix as specified by
    !                  input permutation (used when an ordering has already
    !                  been determined and only the nonzero entries in the
    !                  upper triangle of m are being stored)
    !             4  same as 2 but put diagonal entries at start of each row
    !             5  same as 3 but put diagonal entries at start of each row
    !
    !    flag - integer error flag.  values and their meanings are -
    !               0    no errors detected
    !              9n+k  insufficient storage in md
    !             10n+1  insufficient storage in odrv
    !             11n+1  illegal path specification
    !
    !
    !  conversion from real to double precision
    !
    !    change the real declarations in odrv and sro to double precision
    !    declarations.
    !
    !-----------------------------------------------------------------------
    !
    implicit none

    ! -- dummy variables
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: nja
    integer(I4B), intent(in) :: nsp
    integer(I4B), dimension(n + 1), intent(in) :: ia
    integer(I4B), dimension(nja), intent(in) :: ja
    integer(I4B), dimension(n), intent(inout) :: p
    integer(I4B), dimension(n), intent(inout) :: ip
    integer(I4B), dimension(nsp), intent(inout) :: isp
    integer(I4B), intent(inout) :: flag

    ! -- local
    integer(I4B) :: v
    integer(I4B) :: l
    integer(I4B) :: head
    integer(I4B) :: mmax
    integer(I4B) :: next
    integer(I4B) :: path
    !
    ! set path for finding ordering only
    !
    path = 1
    !
    !
    ! initialize error flag and validate path specification
    flag = 0
    if (path < 1 .or. 5 < path) go to 111
    !
    ! find minimum degree ordering
    mmax = (nsp - n) / 2
    v = 1
    l = v + mmax
    head = l + mmax
    next = head + n
    if (mmax < n) go to 110
    !
    call ims_md(n, nja, ia, ja, mmax, isp(v), isp(l), isp(head), p, &
                ip, isp(v), flag)
    if (flag .ne. 0) go to 100
    !
    return
    !
    ! ** error -- error detected in md
    !             flag = 9 * n + vi from routine mdi.
    !
100 return
    ! ** error -- insufficient storage
110 flag = 10 * n + 1
    return
    ! ** error -- illegal path specified
111 flag = 11 * n + 1
    return
  end subroutine ims_odrv

  subroutine ims_md(n, nja, ia, ja, mmax, v, l, head, last, next, &
                    mark, flag)
    !
    !*****************************************************************
    !  ims_md -- minimum degree algorithm (based on element model)
    !*****************************************************************
    !
    !  description
    !
    !    ims_md finds a minimum degree ordering of the rows and
    !    columns of a general sparse matrix m stored in (ia,ja,a)
    !    format. when the structure of m is nonsymmetric, the ordering
    !    is that obtained for the symmetric matrix  m + m-transpose.
    !
    !
    !  additional parameters
    !
    !    mmax  - declared dimension of the one-dimensional arrays v and l.
    !           mmax must be at least  n+2k,  where k is the number of
    !           nonzeroes in the strict upper triangle of m
    !
    !    v    - integer one-dimensional work array.  dimension = mmax
    !
    !    l    - integer one-dimensional work array.  dimension = mmax
    !
    !    head - integer one-dimensional work array.  dimension = n
    !
    !    last - integer one-dimensional array used to return the permutation
    !           of the rows and columns of m corresponding to the minimum
    !           degree ordering.  dimension = n
    !
    !    next - integer one-dimensional array used to return the inverse of
    !           the permutation returned in last.  dimension = n
    !
    !    mark - integer one-dimensional work array (may be the same as v).
    !           dimension = n
    !
    !    flag - integer error flag.  values and their meanings are -
    !             0      no errors detected
    !             11n+1  insufficient storage in md
    !
    !
    !  definitions of internal parameters
    !
    !    ---------+---------------------------------------------------------
    !    v(s)     - value field of list entry
    !    ---------+---------------------------------------------------------
    !    l(s)     - link field of list entry  (0 =) end of list)
    !    ---------+---------------------------------------------------------
    !    l(vi)    - pointer to element list of uneliminated vertex vi
    !    ---------+---------------------------------------------------------
    !    l(ej)    - pointer to boundary list of active element ej
    !    ---------+---------------------------------------------------------
    !    head(d)  - vj =) vj head of d-list d
    !             -  0 =) no vertex in d-list d
    !
    !
    !             -                  vi uneliminated vertex
    !             -          vi in ek           -       vi not in ek
    !    ---------+-----------------------------+---------------------------
    !    next(vi) - undefined but nonnegative   - vj =) vj next in d-list
    !             -                             -  0 =) vi tail of d-list
    !    ---------+-----------------------------+---------------------------
    !    last(vi) - (not set until mdp)         - -d =) vi head of d-list d
    !             --vk =) compute degree        - vj =) vj last in d-list
    !             - ej =) vi prototype of ej    -  0 =) vi not in any d-list
    !             -  0 =) do not compute degree -
    !    ---------+-----------------------------+---------------------------
    !    mark(vi) - mark(vk)                    - nonneg. tag .lt. mark(vk)
    !
    !
    !             -                   vi eliminated vertex
    !             -      ei active element      -           otherwise
    !    ---------+-----------------------------+---------------------------
    !    next(vi) - -j =) vi was j-th vertex    - -j =) vi was j-th vertex
    !             -       to be eliminated      -       to be eliminated
    !    ---------+-----------------------------+---------------------------
    !    last(vi) -  m =) size of ei = m        - undefined
    !    ---------+-----------------------------+---------------------------
    !    mark(vi) - -m =) overlap count of ei   - undefined
    !             -       with ek = m           -
    !             - otherwise nonnegative tag   -
    !             -       .lt. mark(vk)         -
    !
    !-----------------------------------------------------------------------
    !
    implicit none

    ! -- dummy variables
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: nja
    integer(I4B), dimension(n + 1), intent(in) :: ia
    integer(I4B), dimension(nja), intent(in) :: ja
    integer(I4B), intent(in) :: mmax
    integer(I4B), dimension(mmax), intent(inout) :: v
    integer(I4B), dimension(mmax), intent(inout) :: l
    integer(I4B), dimension(n), intent(inout) :: head
    integer(I4B), dimension(n), intent(inout) :: last
    integer(I4B), dimension(n), intent(inout) :: next
    integer(I4B), dimension(n), intent(inout) :: mark
    integer(I4B), intent(inout) :: flag

    ! -- local
    integer(I4B) :: tag
    integer(I4B) :: dmin
    integer(I4B), pointer :: vk
    integer(I4B), pointer :: ek
    integer(I4B) :: tail
    integer(I4B) :: k

    allocate (vk)
    ek => vk
    !
    ! initialization
    tag = 0
    call ims_mdi(n, nja, ia, ja, mmax, v, l, head, last, next, &
                 mark, tag, flag)
    if (flag .ne. 0) return
    !
    k = 0
    dmin = 1
    !
    ! while  k .lt. n  do
1   if (k >= n) go to 4
    !
    ! search for vertex of minimum degree
2   if (head(dmin) > 0) go to 3
    dmin = dmin + 1
    go to 2
    !
    ! remove vertex vk of minimum degree from degree list
3   vk = head(dmin)
    head(dmin) = next(vk)
    if (head(dmin) > 0) last(head(dmin)) = -dmin
    !
    ! number vertex vk, adjust tag, and tag vk
    k = k + 1
    next(vk) = -k
    last(ek) = dmin - 1
    tag = tag + last(ek)
    mark(vk) = tag
    !
    ! form element ek from uneliminated neighbors of vk
    call ims_mdm(n, mmax, vk, tail, v, l, last, next, mark)
    !
    ! purge inactive elements and do mass elimination
    call ims_mdp(n, mmax, k, ek, tail, v, l, head, last, next, mark)
    !
    ! update degrees of uneliminated vertices in ek
    call ims_mdu(n, mmax, ek, dmin, v, l, head, last, next, mark)
    !
    go to 1
    !
    ! generate inverse permutation from permutation
4   do k = 1, n
      next(k) = -next(k)
      last(next(k)) = k
    end do
    !
    deallocate (vk)
    !
    return
  end subroutine ims_md

  subroutine ims_mdi(n, nja, ia, ja, mmax, v, l, head, last, next, &
                     mark, tag, flag)
    !
    !***********************************************************************
    !  ims_mdi -- initialization
    !***********************************************************************
    implicit none

    ! -- dummy variables
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: nja
    integer(I4B), dimension(n + 1), intent(in) :: ia
    integer(I4B), dimension(nja), intent(in) :: ja
    integer(I4B), intent(in) :: mmax
    integer(I4B), dimension(mmax), intent(inout) :: v
    integer(I4B), dimension(mmax), intent(inout) :: l
    integer(I4B), dimension(n), intent(inout) :: head
    integer(I4B), dimension(n), intent(inout) :: last
    integer(I4B), dimension(n), intent(inout) :: next
    integer(I4B), dimension(n), intent(inout) :: mark
    integer(I4B), intent(in) :: tag
    integer(I4B), intent(inout) :: flag

    ! -- local
    integer(I4B) :: sfs
    integer(I4B) :: vi
    integer(I4B) :: dvi
    integer(I4B) :: vj
    integer(I4B) :: jmin
    integer(I4B) :: jmax
    integer(I4B) :: j
    integer(I4B) :: lvk
    integer(I4B) :: kmax
    integer(I4B) :: k
    integer(I4B) :: nextvi
    integer(I4B) :: ieval
    !
    ! initialize degrees, element lists, and degree lists
    do vi = 1, n
      mark(vi) = 1
      l(vi) = 0
      head(vi) = 0
    end do
    sfs = n + 1
    !
    ! create nonzero structure
    ! for each nonzero entry a(vi,vj)
    louter: do vi = 1, n
      jmin = ia(vi)
      jmax = ia(vi + 1) - 1
      if (jmin > jmax) cycle louter
      linner1: do j = jmin, jmax !5
        vj = ja(j)
        !if (vj-vi) 2, 5, 4
        ieval = vj - vi
        if (ieval == 0) cycle linner1 !5
        if (ieval > 0) go to 4
        !
        ! if a(vi,vj) is in strict lower triangle
        ! check for previous occurrence of a(vj,vi)
        lvk = vi
        kmax = mark(vi) - 1
        if (kmax == 0) go to 4
        linner2: do k = 1, kmax
          lvk = l(lvk)
          if (v(lvk) == vj) cycle linner1 !5
        end do linner2
        ! for unentered entries a(vi,vj)
4       if (sfs >= mmax) go to 101
        !
        ! enter vj in element list for vi
        mark(vi) = mark(vi) + 1
        v(sfs) = vj
        l(sfs) = l(vi)
        l(vi) = sfs
        sfs = sfs + 1
        !
        ! enter vi in element list for vj
        mark(vj) = mark(vj) + 1
        v(sfs) = vi
        l(sfs) = l(vj)
        l(vj) = sfs
        sfs = sfs + 1
      end do linner1
    end do louter
    !
    ! create degree lists and initialize mark vector
    do vi = 1, n
      dvi = mark(vi)
      next(vi) = head(dvi)
      head(dvi) = vi
      last(vi) = -dvi
      nextvi = next(vi)
      if (nextvi > 0) last(nextvi) = vi
      mark(vi) = tag
    end do
    !
    return
    !
    ! ** error-  insufficient storage
101 flag = 9 * n + vi
    return
  end subroutine ims_mdi

  subroutine ims_mdm(n, mmax, vk, tail, v, l, last, next, mark)
    !
    !***********************************************************************
    !  ims_mdm -- form element from uneliminated neighbors of vk
    !***********************************************************************
    implicit none

    ! -- dummy variables
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: mmax
    integer(I4B), intent(in) :: vk
    integer(I4B), intent(inout) :: tail
    integer(I4B), dimension(mmax), intent(inout) :: v
    integer(I4B), dimension(mmax), intent(inout) :: l
    integer(I4B), dimension(n), intent(inout) :: last
    integer(I4B), dimension(n), intent(inout) :: next
    integer(I4B), dimension(n), intent(inout) :: mark

    ! -- local
    integer(I4B) :: tag
    integer(I4B) :: s
    integer(I4B) :: ls
    integer(I4B), pointer :: vs
    integer(I4B), pointer :: es
    integer(I4B) :: b
    integer(I4B) :: lb
    integer(I4B) :: vb
    integer(I4B) :: blp
    integer(I4B) :: blpmax

    allocate (vs)
    es => vs
    !
    ! initialize tag and list of uneliminated neighbors
    tag = mark(vk)
    tail = vk
    !
    ! for each vertex/element vs/es in element list of vk
    ls = l(vk)
1   s = ls
    if (s == 0) go to 5
    ls = l(s)
    vs = v(s)
    if (next(vs) < 0) go to 2
    !
    ! if vs is uneliminated vertex, then tag and append to list of
    ! uneliminated neighbors
    mark(vs) = tag
    l(tail) = s
    tail = s
    go to 4
    !
    ! if es is active element, then ...
    ! for each vertex vb in boundary list of element es
2   lb = l(es)
    blpmax = last(es)
    louter: do blp = 1, blpmax !3
      b = lb
      lb = l(b)
      vb = v(b)
      !
      ! if vb is untagged vertex, then tag and append to list of
      ! uneliminated neighbors
      if (mark(vb) >= tag) cycle louter !3
      mark(vb) = tag
      l(tail) = b
      tail = b
    end do louter
    !
    ! mark es inactive
    mark(es) = tag
    !
4   go to 1
    !
    ! terminate list of uneliminated neighbors
5   l(tail) = 0
    !
    deallocate (vs)
    !
    return
  end subroutine ims_mdm

  subroutine ims_mdp(n, mmax, k, ek, tail, v, l, head, last, next, mark)
    !
    !***********************************************************************
    !  ims_mdp -- purge inactive elements and do mass elimination
    !***********************************************************************
    implicit none

    ! -- dummy variables
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: mmax
    integer(I4B), intent(inout) :: k
    integer(I4B), intent(in) :: ek
    integer(I4B), intent(inout) :: tail
    integer(I4B), dimension(mmax), intent(inout) :: v
    integer(I4B), dimension(mmax), intent(inout) :: l
    integer(I4B), dimension(n), intent(inout) :: head
    integer(I4B), dimension(n), intent(inout) :: last
    integer(I4B), dimension(n), intent(inout) :: next
    integer(I4B), dimension(n), intent(inout) :: mark

    ! -- local
    integer(I4B) :: tag
    integer(I4B) :: free
    integer(I4B) :: li
    integer(I4B) :: vi
    integer(I4B) :: lvi
    integer(I4B) :: evi
    integer(I4B) :: s
    integer(I4B) :: ls
    integer(I4B) :: es
    integer(I4B) :: ilp
    integer(I4B) :: ilpmax
    integer(I4B) :: i
    !
    ! initialize tag
    tag = mark(ek)
    !
    ! for each vertex vi in ek
    li = ek
    ilpmax = last(ek)
    if (ilpmax <= 0) go to 12
    louter: do ilp = 1, ilpmax !11
      i = li
      li = l(i)
      vi = v(li)
      !
      ! remove vi from degree list
      if (last(vi) == 0) go to 3
      if (last(vi) > 0) go to 1
      head(-last(vi)) = next(vi)
      go to 2
1     next(last(vi)) = next(vi)
2     if (next(vi) > 0) last(next(vi)) = last(vi)
      !
      ! remove inactive items from element list of vi
3     ls = vi
4     s = ls
      ls = l(s)
      if (ls == 0) go to 6
      es = v(ls)
      if (mark(es) < tag) go to 5
      free = ls
      l(s) = l(ls)
      ls = s
5     go to 4
      !
      ! if vi is interior vertex, then remove from list and eliminate

6     lvi = l(vi)
      if (lvi .ne. 0) go to 7
      l(i) = l(li)
      li = i
      !
      k = k + 1
      next(vi) = -k
      last(ek) = last(ek) - 1
      cycle louter !11
      !
      ! else ...
      ! classify vertex vi
7     if (l(lvi) .ne. 0) go to 9
      evi = v(lvi)
      if (next(evi) >= 0) go to 9
      if (mark(evi) < 0) go to 8
      !
      ! if vi is prototype vertex, then mark as such, initialize
      ! overlap count for corresponding element, and move vi to end
      ! of boundary list
      last(vi) = evi
      mark(evi) = -1
      l(tail) = li
      tail = li
      l(i) = l(li)
      li = i
      go to 10
      !
      ! else if vi is duplicate vertex, then mark as such and adjust
      ! overlap count for corresponding element
8     last(vi) = 0
      mark(evi) = mark(evi) - 1
      go to 10
      !
      ! else mark vi to compute degree
9     last(vi) = -ek
      !
      ! insert ek in element list of vi
10    v(free) = ek
      l(free) = l(vi)
      l(vi) = free
    end do louter !11
    !
    ! terminate boundary list
12  l(tail) = 0
    !
    return
  end subroutine ims_mdp

  subroutine ims_mdu(n, mmax, ek, dmin, v, l, head, last, next, mark)
    !
    !***********************************************************************
    !  ims_mdu -- update degrees of uneliminated vertices in ek
    !***********************************************************************
    implicit none

    ! -- dummy variables
    integer(I4B), intent(in) :: n
    integer(I4B), intent(in) :: mmax
    integer(I4B), intent(in) :: ek
    integer(I4B), intent(inout) :: dmin
    integer(I4B), dimension(mmax), intent(inout) :: v
    integer(I4B), dimension(mmax), intent(inout) :: l
    integer(I4B), dimension(n), intent(inout) :: head
    integer(I4B), dimension(n), intent(inout) :: last
    integer(I4B), dimension(n), intent(inout) :: next
    integer(I4B), dimension(n), intent(inout) :: mark

    ! -- local
    integer(I4B) :: tag
    integer(I4B) :: vi
    integer(I4B) :: evi
    integer(I4B) :: dvi
    integer(I4B) :: s
    integer(I4B), pointer :: vs
    integer(I4B), pointer :: es
    integer(I4B) :: b
    integer(I4B) :: vb
    integer(I4B) :: ilp
    integer(I4B) :: ilpmax
    integer(I4B) :: blp
    integer(I4B) :: blpmax
    integer(I4B) :: i

    allocate (vs)
    es => vs
    !
    ! initialize tag
    tag = mark(ek) - last(ek)
    !
    ! for each vertex vi in ek
    i = ek
    ilpmax = last(ek)
    if (ilpmax <= 0) go to 11
    louter: do ilp = 1, ilpmax !10
      i = l(i)
      vi = v(i)
      !if (last(vi))  1, 10, 8
      if (last(vi) == 0) cycle louter !10
      if (last(vi) > 0) goto 8
      !
      ! if vi neither prototype nor duplicate vertex, then merge elements
      ! to compute degree
      tag = tag + 1
      dvi = last(ek)
      !
      ! for each vertex/element vs/es in element list of vi
      s = l(vi)
2     s = l(s)
      if (s == 0) go to 9
      vs = v(s)
      if (next(vs) < 0) go to 3
      !
      ! if vs is uneliminated vertex, then tag and adjust degree
      mark(vs) = tag
      dvi = dvi + 1
      go to 5
      !
      ! if es is active element, then expand
      ! check for outmatched vertex
3     if (mark(es) < 0) go to 6
      !
      ! for each vertex vb in es
      b = es
      blpmax = last(es)
      linner: do blp = 1, blpmax !4
        b = l(b)
        vb = v(b)
        !
        ! if vb is untagged, then tag and adjust degree
        if (mark(vb) >= tag) cycle linner !4
        mark(vb) = tag
        dvi = dvi + 1
      end do linner !4
      !
5     go to 2
      !
      ! else if vi is outmatched vertex, then adjust overlaps but do not
      ! compute degree
6     last(vi) = 0
      mark(es) = mark(es) - 1
7     s = l(s)
      if (s == 0) cycle louter !10
      es = v(s)
      if (mark(es) < 0) mark(es) = mark(es) - 1
      go to 7
      !
      ! else if vi is prototype vertex, then calculate degree by
      ! inclusion/exclusion and reset overlap count
8     evi = last(vi)
      dvi = last(ek) + last(evi) + mark(evi)
      mark(evi) = 0
      !
      ! insert vi in appropriate degree list
9     next(vi) = head(dvi)
      head(dvi) = vi
      last(vi) = -dvi
      if (next(vi) > 0) last(next(vi)) = vi
      if (dvi < dmin) dmin = dvi
      !
    end do louter !10
    !
11  deallocate (vs)
    !
    return
  end subroutine ims_mdu

end module IMSReorderingModule
