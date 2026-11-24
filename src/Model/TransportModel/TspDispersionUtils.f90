!> @brief Dispersion coefficients shared by solute dispersion and
!! heat conduction
!<
module TspDispersionUtilsModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DHALF, DONE
  use BaseDisModule, only: DisBaseType

  implicit none
  private
  public :: calc_dispcoef

contains

  !> @brief Calculate the dispersion coefficient for each connection
  !!
  !! The coefficient combines the half-cell conductances of cells n
  !! and m, each using the effective dispersion normal to the shared
  !! face and the saturated flow area. Connections to inactive or
  !! unsaturated cells get a coefficient of zero.
  !<
  subroutine calc_dispcoef(dis, ibdsat, sat, d11, d22, d33, &
                           angle1, angle2, angle3, dispcoef)
    ! modules
    use HGeoUtilModule, only: hyeff
    ! dummy
    class(DisBaseType) :: dis !< discretization
    integer(I4B), dimension(:), intent(in) :: ibdsat !< nonzero if cell has saturation
    real(DP), dimension(:), intent(in) :: sat !< cell saturation
    real(DP), dimension(:), intent(in) :: d11 !< dispersion coefficient
    real(DP), dimension(:), intent(in) :: d22 !< dispersion coefficient
    real(DP), dimension(:), intent(in) :: d33 !< dispersion coefficient
    real(DP), dimension(:), intent(in) :: angle1 !< rotation angle 1
    real(DP), dimension(:), intent(in) :: angle2 !< rotation angle 2
    real(DP), dimension(:), intent(in) :: angle3 !< rotation angle 3
    real(DP), dimension(:), intent(inout) :: dispcoef !< (size:njas) dispersion coefficient
    ! local
    integer(I4B) :: nodes, n, m, ipos
    real(DP) :: clnm, clmn, dn, dm
    real(DP) :: vg1, vg2, vg3
    integer(I4B) :: ihc, isympos
    integer(I4B) :: iavgmeth
    real(DP) :: satn, satm
    real(DP) :: cond, cn, cm, denom
    real(DP) :: anm, amn

    ! use arithmetic averaging for effective dispersion
    iavgmeth = 1

    nodes = size(d11)
    do n = 1, nodes
      if (ibdsat(n) == 0) cycle
      do ipos = dis%con%ia(n) + 1, dis%con%ia(n + 1) - 1
        m = dis%con%ja(ipos)
        if (m < n) cycle
        isympos = dis%con%jas(ipos)
        dispcoef(isympos) = DZERO
        if (ibdsat(m) == 0) cycle

        clnm = dis%con%cl1(isympos)
        clmn = dis%con%cl2(isympos)
        ihc = dis%con%ihc(isympos)
        satn = sat(n)
        satm = sat(m)

        ! effective dispersion for n and m normal to the shared face
        call dis%connection_normal(n, m, ihc, vg1, vg2, vg3, ipos)
        dn = hyeff(d11(n), d22(n), d33(n), &
                   angle1(n), angle2(n), angle3(n), &
                   vg1, vg2, vg3, iavgmeth)
        dm = hyeff(d11(m), d22(m), d33(m), &
                   angle1(m), angle2(m), angle3(m), &
                   vg1, vg2, vg3, iavgmeth)

        call dis%get_flow_area(n, m, ipos, satn, satm, anm, amn)
        if (ihc == 0) then
          clnm = satn * (dis%top(n) - dis%bot(n)) * DHALF
          clmn = satm * (dis%top(m) - dis%bot(m)) * DHALF

          ! no vertical dispersion if either cell is dry, or if the
          ! higher-numbered cell is partially saturated
          if (satn == DZERO) then
            anm = DZERO
          else if (n > m .and. satn < DONE) then
            anm = DZERO
          end if
          if (satm == DZERO) then
            anm = DZERO
          else if (m > n .and. satm < DONE) then
            anm = DZERO
          end if
          amn = anm
        else
          ! no horizontal dispersion to or from a dry cell
          if (satn == DZERO .or. satm == DZERO) then
            anm = DZERO
            amn = DZERO
          end if
        end if

        ! combine the two half-cell conductances
        cn = DZERO
        if (clnm > DZERO) cn = dn * anm / clnm
        cm = DZERO
        if (clmn > DZERO) cm = dm * amn / clmn
        denom = cn + cm
        if (denom > DZERO) then
          cond = cn * cm / denom
        else
          cond = DZERO
        end if
        dispcoef(isympos) = cond
      end do
    end do
  end subroutine calc_dispcoef

end module TspDispersionUtilsModule
