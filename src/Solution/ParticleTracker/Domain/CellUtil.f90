module CellUtilModule
  use KindModule, only: I4B, DP
  implicit none

  private
  public :: cell_poly_to_rect
  public :: cell_poly_to_quad

contains

  !> @brief Convert CellPoly representation to CellRect.
  !! Assumes the conversion is possible.
  subroutine cell_poly_to_rect(poly, rect)
    use ConstantsModule, only: DONE
    use CellRectModule, only: CellRectType, create_cell_rect
    use CellPolyModule, only: CellPolyType
    use CellDefnModule, only: CellDefnType
    ! dummy
    type(CellPolyType), intent(in), pointer :: poly
    type(CellRectType), intent(inout), pointer :: rect
    ! local
    type(CellDefnType), pointer :: defn
    integer(I4B) :: ipv, ipv1, ipv2, ipv3, ipv4
    integer(I4B), dimension(4) :: ipvnxt = (/2, 3, 4, 1/)
    real(DP) :: x1, y1, x2, y2, x4, y4
    real(DP) :: dx2, dy2, dx4, dy4, areax, areay, areaz
    real(DP) :: xOrigin, yOrigin, zOrigin, dx, dy, dz, sinrot, cosrot
    real(DP) :: factor, term

    call create_cell_rect(rect)
    defn => poly%defn
    ! Translate and rotate the rectangular cell into local coordinates
    ! with x varying from 0 to dx and y varying from 0 to dy. Choose the
    ! "south-west" vertex to be the local origin so that the rotation
    ! angle is zero if the cell already aligns with the model x and y
    ! coordinates. The "southwest" vertex is found by finding the vertex
    ! formed either by (1) an edge over which x decreases (going
    ! clockwise) followed by an edge over which x does not increase, or
    ! by (2) an edge over which y does not decrease (again going
    ! clockwise) followed by an edge over which y increases. In the end,
    ! ipv1 is the local vertex number (within the cell, taking a value
    ! of 1, 2, 3, or 4) of the southwest vertex, and ipv2, ipv3, and
    ! ipv4 are the local vertex numbers of the remaining three vertices
    ! going clockwise.
    do ipv = 1, 4
      ipv4 = ipv
      ipv1 = ipvnxt(ipv4)
      x4 = defn%polyvert(1, ipv4)
      y4 = defn%polyvert(2, ipv4)
      x1 = defn%polyvert(1, ipv1)
      y1 = defn%polyvert(2, ipv1)
      if (x1 .lt. x4) then
        ipv2 = ipvnxt(ipv1)
        x2 = defn%polyvert(1, ipv2)
        if (x2 .le. x1) then
          y2 = defn%polyvert(2, ipv2)
          exit
        end if
      else if (y1 .ge. y4) then
        ipv2 = ipvnxt(ipv1)
        y2 = defn%polyvert(2, ipv2)
        if (y2 .gt. y1) then
          x2 = defn%polyvert(1, ipv2)
          exit
        end if
      end if
    end do
    ipv3 = ipvnxt(ipv2)

    ! Compute upper bounds on the local coordinates (the rectangular
    ! dimensions of the cell) and the sine and cosine of the rotation
    ! angle, and store local origin information
    xOrigin = x1
    yOrigin = y1
    zOrigin = defn%bot
    dx2 = x2 - xOrigin
    dy2 = y2 - yOrigin
    dx4 = x4 - xOrigin
    dy4 = y4 - yOrigin
    dx = dsqrt(dx4 * dx4 + dy4 * dy4)
    dy = dsqrt(dx2 * dx2 + dy2 * dy2)
    dz = defn%top - zOrigin
    sinrot = dy4 / dx
    cosrot = dx4 / dx
    rect%defn = poly%defn
    rect%dx = dx
    rect%dy = dy
    rect%dz = dz
    rect%sinrot = sinrot
    rect%cosrot = cosrot
    rect%xOrigin = xOrigin
    rect%yOrigin = yOrigin
    rect%zOrigin = zOrigin
    rect%ipvOrigin = ipv1

    ! Compute (unscaled) cell edge velocities from face flows
    areax = dx * dz
    areay = dy * dz
    areaz = dx * dy
    factor = DONE / (defn%retfactor * defn%porosity)
    term = factor / areax
    rect%vx1 = defn%faceflow(ipv1) * term
    rect%vx2 = -defn%faceflow(ipv3) * term
    term = factor / areay
    rect%vy1 = defn%faceflow(ipv4) * term
    rect%vy2 = -defn%faceflow(ipv2) * term
    term = factor / areaz
    rect%vz1 = defn%faceflow(6) * term
    rect%vz2 = -defn%faceflow(7) * term
  end subroutine cell_poly_to_rect

  !> @brief Convert CellPoly representation to CellRectQuad.
  !! Assumes the conversion is possible.
  subroutine cell_poly_to_quad(poly, quad)
    use CellRectQuadModule, only: CellRectQuadType, create_cell_rect_quad
    use CellPolyModule, only: CellPolyType
    use MathUtilModule, only: mod_offset
    ! dummy
    type(CellPolyType), intent(in), pointer :: poly
    type(CellRectQuadType), intent(inout), pointer :: quad
    ! local
    integer(I4B) :: i, irv, isc
    real(DP) :: qhalf, qdisttopbot, q1, q2, q4

    call create_cell_rect_quad(quad)
    call quad%init_from(poly%defn)
    ! Translate and rotate the rect-quad cell into local coordinates with
    ! x varying from 0 to dx and y varying from 0 to dy. Choose the "south-
    ! west" rectangle vertex to be the local origin so that the rotation
    ! angle is zero if the cell already aligns with the model x and y
    ! coordinates.
    quad%irvOrigin = quad%get_rect_ivert_sw()
    call quad%get_rect_dim_rot()

    ! Set the external and internal face flows used for subcells
    do i = 0, 3
      irv = mod_offset(i + quad%irvOrigin, 4, 1)
      isc = mod_offset(i + 3, 4, 1)
      if (.not. quad%face_is_refined(irv)) then
        qhalf = 5d-1 * quad%get_rect_flow(1, irv)
        quad%qextl2(isc) = qhalf
        isc = mod_offset(isc + 1, 4, 1)
        quad%qextl1(isc) = qhalf
      else
        quad%qextl2(isc) = quad%get_rect_flow(1, irv)
        isc = mod_offset(isc + 1, 4, 1)
        quad%qextl1(isc) = quad%get_rect_flow(2, irv)
      end if
    end do
    qdisttopbot = 2.5d-1 * (quad%defn%get_distflow() &
                            + quad%defn%get_botflow() &
                            + quad%defn%get_topflow())
    q1 = qdisttopbot + quad%qextl1(1) + quad%qextl2(1)
    q2 = qdisttopbot + quad%qextl1(2) + quad%qextl2(2)
    q4 = qdisttopbot + quad%qextl1(4) + quad%qextl2(4)
    quad%qintl(1) = -5d-1 * (q1 + 5d-1 * (q2 - q4))
    quad%qintl(2) = quad%qintl(1) + q1
    quad%qintl(3) = quad%qintl(2) + q2
    quad%qintl(4) = quad%qintl(1) - q4
    quad%qintl(5) = quad%qintl(1)
  end subroutine cell_poly_to_quad

end module CellUtilModule
