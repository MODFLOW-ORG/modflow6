module GeomUtilModule
  use KindModule, only: I4B, DP, LGP
  use ErrorUtilModule, only: pstop
  use ConstantsModule, only: DZERO, DONE, DTWO, DHALF, &
                             DONETHIRD, DEP3

  implicit none
  private
  public :: between, point_in_polygon, &
            get_node, get_ijk, get_jk, &
            skew, transform, compose, &
            area, polygon_extent, shared_face, clamp_bary
contains

  !> @brief Check if a value is between two other values (inclusive).
  logical function between(x, a, b)
    real(DP), intent(in) :: x, a, b
    between = ((x >= a .and. x <= b) .or. (x <= a .and. x >= b))
  end function between

  !> @brief Check if a point is within a polygon.
  !!
  !! Vertices and edge points are considered in the polygon. By
  !! default, a point must lie exactly on an edge to be considered
  !! within the polygon. A tolerance may be specified to instead
  !! accept coordinates within tol of an edge.
  !!
  !! Adapted from https://stackoverflow.com/a/63436180/6514033,
  !<
  logical function point_in_polygon(x, y, poly, tol)
    ! dummy
    real(DP), intent(in) :: x !< x point coordinate
    real(DP), intent(in) :: y !< y point coordinate
    real(DP), allocatable, intent(in) :: poly(:, :) !< polygon vertices (column-major indexing)
    real(DP), intent(in), optional :: tol !< tolerance (default 0)
    ! local
    integer(I4B) :: i, ii, num_verts
    real(DP) :: xa, xb, ya, yb, c
    real(DP) :: ltol

    ltol = DZERO
    if (present(tol)) ltol = tol

    point_in_polygon = .false.
    num_verts = size(poly, 2)
    xa = poly(1, num_verts)
    ya = poly(2, num_verts)

    do i = 0, num_verts - 1
      ii = mod(i, num_verts) + 1
      xb = poly(1, ii)
      yb = poly(2, ii)

      if ((x == xa .and. y == ya) .or. &
          (x == xb .and. y == yb)) then
        ! on vertex
        point_in_polygon = .true.
        exit
      else if (ya == yb .and. &
               abs(y - ya) * abs(xb - xa) <= ltol .and. &
               between(x, xa, xb)) then
        ! on (or within tol of) horizontal edge
        point_in_polygon = .true.
        exit
      else if (between(y, ya, yb)) then
        if ((y == ya .and. yb >= ya) .or. &
            (y == yb .and. ya >= yb)) then
          xa = xb
          ya = yb
          cycle
        end if
        ! cross product
        c = (xa - x) * (yb - y) - (xb - x) * (ya - y)
        if (abs(c) <= ltol) then
          ! on (or within tol of) edge
          point_in_polygon = .true.
          exit
        else if ((ya < yb) .eqv. (c > 0)) then
          ! ray intersection
          point_in_polygon = .not. point_in_polygon
        end if
      end if

      xa = xb
      ya = yb
    end do
  end function point_in_polygon

  !> @brief Get node number, given layer, row, and column indices
  !! for a structured grid. If any argument is invalid return -1.
  function get_node(ilay, irow, icol, nlay, nrow, ncol)
    integer(I4B), intent(in) :: ilay, irow, icol, nlay, nrow, ncol
    integer(I4B) :: get_node

    if (nlay > 0 .and. nrow > 0 .and. ncol > 0 .and. &
        ilay > 0 .and. ilay <= nlay .and. &
        irow > 0 .and. irow <= nrow .and. &
        icol > 0 .and. icol <= ncol) then
      get_node = &
        icol + ncol * (irow - 1) + (ilay - 1) * nrow * ncol
    else
      get_node = -1
    end if
  end function get_node

  !> @brief Get row, column and layer indices from node number and grid
  !! dimensions. If nodenumber is invalid, irow, icol, and ilay are -1.
  subroutine get_ijk(nodenumber, nrow, ncol, nlay, irow, icol, ilay)
    ! -- dummy variables
    integer(I4B), intent(in) :: nodenumber
    integer(I4B), intent(in) :: nrow
    integer(I4B), intent(in) :: ncol
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(out) :: irow
    integer(I4B), intent(out) :: icol
    integer(I4B), intent(out) :: ilay
    ! -- local variables
    integer(I4B) :: nodes
    integer(I4B) :: ij

    nodes = nlay * nrow * ncol
    if (nodenumber < 1 .or. nodenumber > nodes) then
      irow = -1
      icol = -1
      ilay = -1
    else
      ilay = (nodenumber - 1) / (ncol * nrow) + 1
      ij = nodenumber - (ilay - 1) * ncol * nrow
      irow = (ij - 1) / ncol + 1
      icol = ij - (irow - 1) * ncol
    end if
  end subroutine get_ijk

  !> @brief Get layer index and within-layer node index from node number
  !! and grid dimensions. If nodenumber is invalid, icpl and ilay are -1.
  subroutine get_jk(nodenumber, ncpl, nlay, icpl, ilay)
    ! -- dummy variables
    integer(I4B), intent(in) :: nodenumber
    integer(I4B), intent(in) :: ncpl
    integer(I4B), intent(in) :: nlay
    integer(I4B), intent(out) :: icpl
    integer(I4B), intent(out) :: ilay
    ! -- local variables
    integer(I4B) :: nodes

    nodes = ncpl * nlay
    if (nodenumber < 1 .or. nodenumber > nodes) then
      icpl = -1
      ilay = -1
    else
      ilay = (nodenumber - 1) / ncpl + 1
      icpl = nodenumber - (ilay - 1) * ncpl
    end if
  end subroutine get_jk

  !> @brief Skew a 2D vector along the x-axis.
  pure function skew(v, s, invert) result(res)
    ! -- dummy
    real(DP), intent(in) :: v(2) !< vector
    real(DP), intent(in) :: s(3) !< skew matrix entries (top left, top right, bottom right)
    logical(LGP), intent(in), optional :: invert
    real(DP) :: res(2)
    ! -- local
    logical(LGP) :: linvert
    real(DP) :: sxx, sxy, syy

    ! -- process optional arguments
    if (present(invert)) then
      linvert = invert
    else
      linvert = .false.
    end if

    sxx = s(1)
    sxy = s(2)
    syy = s(3)
    if (.not. linvert) then
      res(1) = sxx * v(1) + sxy * v(2)
      res(2) = syy * v(2)
    else
      res(2) = v(2) / syy
      res(1) = (v(1) - sxy * res(2)) / sxx
    end if
  end function skew

  !> @brief Apply a 3D translation and optional 2D rotation to coordinates.
  subroutine transform(xin, yin, zin, &
                       xout, yout, zout, &
                       xorigin, yorigin, zorigin, &
                       sinrot, cosrot, &
                       invert)
    ! -- dummy
    real(DP) :: xin, yin, zin !< input coordinates
    real(DP) :: xout, yout, zout !< output coordinates
    real(DP), optional :: xorigin, yorigin, zorigin !< origin coordinates
    real(DP), optional :: sinrot, cosrot !< sine and cosine of rotation
    logical(LGP), optional :: invert !< whether to invert
    ! -- local
    logical(LGP) :: ltranslate, lrotate, linvert
    real(DP) :: x, y
    real(DP) :: lxorigin, lyorigin, lzorigin
    real(DP) :: lsinrot, lcosrot

    ! -- Process option arguments and set defaults and flags
    call defaults(lxorigin, lyorigin, lzorigin, &
                  lsinrot, lcosrot, linvert, &
                  ltranslate, lrotate, &
                  xorigin, yorigin, zorigin, &
                  sinrot, cosrot, invert)

    ! -- Apply transformation or its inverse
    if (.not. linvert) then
      ! -- Apply transformation to coordinates
      if (ltranslate) then
        xout = xin - lxorigin
        yout = yin - lyorigin
        zout = zin - lzorigin
      else
        xout = lxorigin
        yout = lyorigin
        zout = lzorigin
      end if
      if (lrotate) then
        x = xout
        y = yout
        xout = x * lcosrot + y * lsinrot
        yout = -x * lsinrot + y * lcosrot
      end if
    else
      ! -- Apply inverse of transformation to coordinates
      if (lrotate) then
        x = xin * lcosrot - yin * lsinrot
        y = xin * lsinrot + yin * lcosrot
      else
        x = xin
        y = yin
      end if
      if (ltranslate) then
        xout = x + lxorigin
        yout = y + lyorigin
        zout = zin + lzorigin
      end if
    end if
  end subroutine transform

  !> @brief Compose an affine transform onto an existing one, or remove (invert) it.
  subroutine compose(xorigin, yorigin, zorigin, &
                     sinrot, cosrot, &
                     xorigin_new, yorigin_new, zorigin_new, &
                     sinrot_new, cosrot_new, &
                     invert)
    ! -- dummy
    real(DP) :: xorigin, yorigin, zorigin !< cumulative transform T, origin (in/out)
    real(DP) :: sinrot, cosrot !< cumulative transform T, rotation as (sin, cos) (in/out)
    real(DP), optional :: xorigin_new, yorigin_new, zorigin_new !< transform A to compose, origin
    real(DP), optional :: sinrot_new, cosrot_new !< transform A to compose, rotation as (sin, cos)
    logical(LGP), optional :: invert !< whether to remove A rather than compose it
    ! -- local
    logical(LGP) :: ltranslate, lrotate, linvert
    real(DP) :: xa, ya, za, sa, ca !< transform A (origin, then sin/cos of phi)
    real(DP) :: x0, y0, z0, s0, c0 !< incoming T (origin, then sin/cos of its angle)

    ! -- Process option arguments and set defaults and flags
    call defaults(xa, ya, za, sa, ca, linvert, &
                  ltranslate, lrotate, &
                  xorigin_new, yorigin_new, zorigin_new, &
                  sinrot_new, cosrot_new, invert)

    ! -- Copy existing transformation into working copy
    x0 = xorigin
    y0 = yorigin
    z0 = zorigin
    s0 = sinrot
    c0 = cosrot

    ! -- Given transformations A = (a, phi) and T = (t, alpha)
    ! --   forward:  alpha' = alpha + phi,   t' = t + R(alpha) a
    ! --   inverse:  alpha  = alpha' - phi,  t  = t' - R(alpha' - phi) a
    if (.not. linvert) then
      if (lrotate) then
        ! -- alpha' = alpha + phi
        sinrot = ca * s0 + sa * c0
        cosrot = ca * c0 - sa * s0
      end if
      if (ltranslate) then
        ! -- t' = t + R(alpha) a, with R(alpha) built from the incoming (s0, c0)
        xorigin = x0 + (c0 * xa - s0 * ya)
        yorigin = y0 + (s0 * xa + c0 * ya)
        zorigin = z0 + za
      end if
    else
      if (lrotate) then
        ! -- alpha = alpha' - phi
        ! -- (sinrot, cosrot) now hold R(alpha' - phi)
        sinrot = ca * s0 - sa * c0
        cosrot = ca * c0 + sa * s0
      end if
      if (ltranslate) then
        ! -- t = t' - R(alpha' - phi) a. For phi = 0 this is t' - R(alpha') a.
        ! -- if A carries no rotation, (sinrot, cosrot) still hold R(alpha').
        xorigin = x0 - (cosrot * xa - sinrot * ya)
        yorigin = y0 - (sinrot * xa + cosrot * ya)
        zorigin = z0 - za
      end if
    end if
  end subroutine compose

  !> @brief Process arguments and set defaults. Internal use only.
  subroutine defaults(xorigin, yorigin, zorigin, &
                      sinrot, cosrot, &
                      invert, translate, rotate, &
                      xorigin_opt, yorigin_opt, zorigin_opt, &
                      sinrot_opt, cosrot_opt, invert_opt)
    ! -- dummy
    real(DP) :: xorigin, yorigin, zorigin
    real(DP) :: sinrot, cosrot
    logical(LGP) :: invert, translate, rotate
    real(DP), optional :: xorigin_opt, yorigin_opt, zorigin_opt
    real(DP), optional :: sinrot_opt, cosrot_opt
    logical(LGP), optional :: invert_opt

    translate = .false.
    xorigin = DZERO
    if (present(xorigin_opt)) then
      xorigin = xorigin_opt
      translate = .true.
    end if
    yorigin = DZERO
    if (present(yorigin_opt)) then
      yorigin = yorigin_opt
      translate = .true.
    end if
    zorigin = DZERO
    if (present(zorigin_opt)) then
      zorigin = zorigin_opt
      translate = .true.
    end if
    rotate = .false.
    sinrot = DZERO
    cosrot = DONE
    if (present(sinrot_opt)) then
      sinrot = sinrot_opt
      if (present(cosrot_opt)) then
        cosrot = cosrot_opt
      else
        ! -- If sinrot_opt is specified but cosrot_opt is not,
        ! -- default to corresponding non-negative cosrot
        cosrot = dsqrt(DONE - sinrot * sinrot)
      end if
      rotate = .true.
    else if (present(cosrot_opt)) then
      cosrot = cosrot_opt
      ! -- cosrot_opt is specified but sinrot_opt is not, so
      ! -- default to corresponding non-negative sinrot
      sinrot = dsqrt(DONE - cosrot * cosrot)
      rotate = .true.
    end if
    invert = .false.
    if (present(invert_opt)) invert = invert_opt
  end subroutine defaults

  !> @brief Calculate polygon area, with vertices given in CW or CCW order.
  function area(xv, yv, cw) result(a)
    ! dummy
    real(DP), dimension(:), intent(in) :: xv
    real(DP), dimension(:), intent(in) :: yv
    logical(LGP), intent(in), optional :: cw
    ! result
    real(DP) :: a
    integer(I4B) :: s

    if (present(cw)) then
      if (cw) then
        s = 1
      else
        s = -1
      end if
    else
      s = 1
    end if

    a = -DHALF * sum(xv(:) * cshift(yv(:), s) - cshift(xv(:), s) * yv(:))

  end function area

  !> @brief Calculate the maximum distance between two polygon vertices.
  pure function polygon_extent(xv, yv) result(e)
    ! dummy
    real(DP), dimension(:), intent(in) :: xv
    real(DP), dimension(:), intent(in) :: yv
    ! result
    real(DP) :: e
    ! local
    integer(I4B) :: i
    integer(I4B) :: j
    real(DP) :: dx
    real(DP) :: dy

    e = DZERO
    do i = 1, size(xv) - 1
      do j = i + 1, size(xv)
        dx = xv(i) - xv(j)
        dy = yv(i) - yv(j)
        e = max(e, dx * dx + dy * dy)
      end do
    end do
    e = sqrt(e)

  end function polygon_extent

  !> @brief Find the lateral face shared by two cells.
  !!
  !! Find the lateral (x-y plane) face shared by the given cells.
  !! The iface return argument will be 0 if they share no such face,
  !! otherwise the index of the shared face in cell 1's vertex array,
  !! where face N connects vertex N to vertex N + 1 going clockwise.
  !!
  !! Note: assumes the cells are convex and share at most 2 vertices
  !! and that both vertex arrays are oriented clockwise.
  !<
  subroutine shared_face(iverts1, iverts2, iface)
    integer(I4B), dimension(:) :: iverts1
    integer(I4B), dimension(:) :: iverts2
    integer(I4B), intent(out) :: iface
    integer(I4B) :: nv1
    integer(I4B) :: nv2
    integer(I4B) :: il1, iil1
    integer(I4B) :: il2, iil2
    logical(LGP) :: found
    logical(LGP) :: wrapped

    iface = 0
    found = .false.
    nv1 = size(iverts1)
    nv2 = size(iverts2)
    wrapped = iverts1(1) == iverts1(nv1)

    ! Find a vertex shared by the cells, then check the adjacent faces.
    ! If the cells share a face, it must be one of these. When looking
    ! forward in the 1st cell's vertices, look backwards in the 2nd's,
    ! and vice versa, since a clockwise face in cell 1 must correspond
    ! to a counter-clockwise face in cell 2.
    outerloop: do il1 = 1, nv1 - 1
      do il2 = 1, nv2 - 1
        if (iverts1(il1) == iverts2(il2)) then

          iil1 = il1 + 1
          if (il2 == 1) then
            iil2 = nv2
            if (wrapped) iil2 = iil2 - 1
          else
            iil2 = il2 - 1
          end if
          if (iverts1(iil1) == iverts2(iil2)) then
            found = .true.
            iface = il1
            exit outerloop
          end if

          iil2 = il2 + 1
          if (il1 == 1) then
            iil1 = nv1
            if (wrapped) iil1 = iil1 - 1
          else
            iil1 = il1 - 1
          end if
          if (iverts1(iil1) == iverts2(iil2)) then
            found = .true.
            iface = iil1
            exit outerloop
          end if

        end if
      end do
      if (found) exit
    end do outerloop
  end subroutine shared_face

  !> @brief Clamp barycentric coordinates to the interior of a triangle,
  !! with optional padding some minimum distance from any face.
  !!
  !! This routine requires 0 <= tol <= 1/3 and 1 = alpha + beta + gamma.
  !<
  subroutine clamp_bary(alpha, beta, gamma, pad)
    ! dummy
    real(DP), intent(inout) :: alpha
    real(DP), intent(inout) :: beta
    real(DP), intent(out) :: gamma
    real(DP), intent(in), optional :: pad
    ! local
    real(DP) :: lolimit
    real(DP) :: hilimit
    real(DP) :: delta
    real(DP) :: lpad

    if (present(pad)) then
      lpad = pad
      if (pad < DZERO .or. pad > DONETHIRD) &
        call pstop(1, "pad must be between 0 and 1/3, inclusive")
    else
      lpad = DZERO
    end if

    gamma = DONE - alpha - beta
    lolimit = lpad
    hilimit = DONE - DTWO * lpad
    ! Check alpha coordinate against lower limit
    if (alpha < lolimit) then
      ! Alpha is too low, so nudge alpha to lower limit; this is a move
      ! parallel to the "alpha axis," which also changes gamma
      alpha = lolimit
      gamma = DONE - alpha - beta
      ! Check beta coordinate against lower limit (which in this
      ! case is equivalent to checking gamma coordinate against
      ! upper limit)
      if (beta < lolimit) then
        ! Beta is too low (gamma is too high), so nudge beta to lower limit;
        ! this is a move parallel to the "beta axis," which also changes gamma
        beta = lolimit
        gamma = hilimit
        ! Check beta coordinate against upper limit (which in this
        ! case is equivalent to checking gamma coordinate against
        ! lower limit)
      else if (beta > hilimit) then
        ! Beta is too high (gamma is too low), so nudge beta to lower limit;
        ! this is a move parallel to the "beta axis," which also changes gamma
        beta = hilimit
        gamma = lolimit
      end if
    end if
    ! Check beta coordinate against lower limit. (If alpha coordinate
    ! was nudged to lower limit, beta and gamma coordinates have also
    ! been adjusted as necessary to place particle within subcell, and
    ! subsequent checks on beta and gamma will evaluate to false, and
    ! no further adjustments will be made.)
    if (beta < lolimit) then
      ! Beta is too low, so nudge beta to lower limit; this is a move
      ! parallel to the "beta axis," which also changes gamma
      beta = lolimit
      gamma = DONE - alpha - beta
      ! Check alpha coordinate against lower limit (which in this
      ! case is equivalent to checking gamma coordinate against
      ! upper limit)
      if (alpha < lolimit) then
        ! Alpha is too low (gamma is too high), so nudge alpha to lower limit;
        ! this is a move parallel to the "alpha axis," which also changes gamma
        alpha = lolimit
        gamma = hilimit
        ! Check alpha coordinate against upper limit (which in this
        ! case is equivalent to checking gamma coordinate against
        ! lower limit)
      else if (alpha > hilimit) then
        ! Alpha is too high (gamma is too low), so nudge alpha to lower limit;
        ! this is a move parallel to the "alpha axis," which also changes gamma
        alpha = hilimit
        gamma = lolimit
      end if
    end if
    ! Check gamma coordinate against lower limit.(If alpha and/or beta
    ! coordinate was nudged to lower limit, gamma coordinate has also
    ! been adjusted as necessary to place particle within subcell, and
    ! subsequent check on gamma will evaluate to false, and no further
    ! adjustment will be made.)
    if (gamma < lolimit) then
      ! Gamma is too low, so nudge gamma to lower limit; this is a move
      ! parallel to the "gamma axis," which also changes alpha and beta
      delta = DHALF * (lolimit - gamma)
      gamma = lpad
      alpha = alpha - delta
      beta = beta - delta
      ! Check beta coordinate against lower limit (which in this
      ! case is equivalent to checking alpha coordinate against
      ! upper limit)
      if (beta < lolimit) then
        ! Beta is too low (alpha is too high), so nudge beta to lower limit;
        ! this is a move parallel to the "gamma axis," which also changes alpha
        beta = lolimit
        alpha = DONE - beta - gamma
        ! Check beta coordinate against upper limit (which in this
        ! case is equivalent to checking gamma coordinate against
        ! lower limit)
      else if (beta > hilimit) then
        ! Beta is too high (alpha is too low), so nudge beta to lower limit;
        ! this is a move parallel to the "gamma axis," which also changes alpha
        beta = hilimit
        alpha = DONE - beta - gamma
      end if
    end if
  end subroutine clamp_bary

end module GeomUtilModule
