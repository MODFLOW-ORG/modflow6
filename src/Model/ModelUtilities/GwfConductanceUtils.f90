!> @brief This module contains stateless conductance functions
!!
!! This module contains the functions to calculate the horizontal
!! and vertical conductance between two cells that are used in the
!! the node property flow (NPF) package. It also contains functions
!! to calculate the wetted cell fraction. This module does not
!! depend on the NPF package.
!!
!<
module GwfConductanceUtilsModule
  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DHALF, DONE, &
                             DLNLOW, DLNHIGH, &
                             C3D_STAGGERED

  implicit none
  private
  public :: hcond
  public :: vcond
  public :: condmean
  public :: thksatnm
  public :: staggered_thkfrac
  public :: CCOND_HMEAN

  !> @brief enumerator that defines the conductance options
  !<
  ENUM, BIND(C)
    ENUMERATOR :: CCOND_HMEAN = 0 !< Harmonic mean
    ENUMERATOR :: CCOND_LMEAN = 1 !< Logarithmic mean
    ENUMERATOR :: CCOND_AMTLMK = 2 !< Arithmetic-mean thickness and logarithmic-mean hydraulic conductivity
    ENUMERATOR :: CCOND_AMTHMK = 3 !< Arithmetic-mean thickness and harmonic-mean hydraulic conductivity
  END ENUM

contains

  !> @brief Horizontal conductance between two cells
  !!
  !! This function uses a weighted transmissivity in the harmonic mean
  !! conductance calculations. This differs from the MODFLOW-NWT and
  !! MODFLOW-USG conductance calculations for the Newton-Raphson formulation
  !! which use a weighted hydraulic conductivity.
  !<
  function hcond(ibdn, ibdm, ictn, ictm, iupstream, ihc, icellavg, &
                 condsat, hn, hm, satn, satm, hkn, hkm, topn, topm, &
                 botn, botm, cln, clm, fawidth) &
    result(condnm)
    ! return variable
    real(DP) :: condnm !< horizontal conductance between two cells
    ! dummy
    integer(I4B), intent(in) :: ibdn !< cell n active flag
    integer(I4B), intent(in) :: ibdm !< cell m active flag
    integer(I4B), intent(in) :: ictn !< cell n convertible cell flag
    integer(I4B), intent(in) :: ictm !< cell m convertible cell flag
    integer(I4B), intent(in) :: iupstream !< flag for upstream weighting
    integer(I4B), intent(in) :: ihc !< connection type
    integer(I4B), intent(in) :: icellavg !< cell averaging option
    real(DP), intent(in) :: condsat !< saturated conductance
    real(DP), intent(in) :: hn !< cell n head
    real(DP), intent(in) :: hm !< cell m head
    real(DP), intent(in) :: satn !< cell n wetted cell fraction
    real(DP), intent(in) :: satm !< cell m wetted cell fraction
    real(DP), intent(in) :: hkn !< horizontal hydraulic conductivity for cell n (in the direction of cell m)
    real(DP), intent(in) :: hkm !< horizontal hydraulic conductivity for cell m (in the direction of cell n)
    real(DP), intent(in) :: topn !< top of cell n
    real(DP), intent(in) :: topm !< top of cell m
    real(DP), intent(in) :: botn !< bottom of cell n
    real(DP), intent(in) :: botm !< bottom of cell m
    real(DP), intent(in) :: cln !< distance from the center of cell n to the shared face with cell m
    real(DP), intent(in) :: clm !< distance from the center of cell m to the shared face with cell n
    real(DP), intent(in) :: fawidth !< width of cell perpendicular to flow

    ! n or m is inactive
    if (ibdn == 0 .or. ibdm == 0) then
      condnm = DZERO
      ! both cells are non-convertible
    elseif (ictn == 0 .and. ictm == 0) then
      condnm = condsat
    else if (iupstream == 1) then
      condnm = convertible_upstream(hn, hm, satn, satm, condsat)
    else
      condnm = convertible_standard(ihc, icellavg, &
                                    satn, satm, hkn, hkm, &
                                    topn, topm, botn, botm, &
                                    cln, clm, fawidth)
    end if
  end function hcond

  !> @brief Convertible cell(s) with upstream weighted horizontal conductance
  !<
  function convertible_upstream(hn, hm, satn, satm, condsat) &
    result(condnm)
    ! return variable
    real(DP) :: condnm !< horizontal conductance between two cells
    ! dummy
    real(DP), intent(in) :: condsat !< saturated conductance
    real(DP), intent(in) :: hn !< cell n head
    real(DP), intent(in) :: hm !< cell m head
    real(DP), intent(in) :: satn !< cell n wetted cell fraction
    real(DP), intent(in) :: satm !< cell m wetted cell fraction
    !  local
    real(DP) :: sat_up

    if (hn > hm) then
      sat_up = satn
    else
      sat_up = satm
    end if
    condnm = sat_up * condsat
  end function convertible_upstream

  !> @brief Convertible cell(s) with standard weighted horizontal conductance
  !<
  function convertible_standard(ihc, icellavg, satn, satm, hkn, hkm, &
                                topn, topm, botn, botm, cln, clm, fawidth) &
    result(condnm)
    ! return variable
    real(DP) :: condnm !< horizontal conductance between two cells
    ! dummy
    integer(I4B), intent(in) :: ihc !< connection type
    integer(I4B), intent(in) :: icellavg !< cell averaging option
    real(DP), intent(in) :: satn !< cell n wetted cell fraction
    real(DP), intent(in) :: satm !< cell m wetted cell fraction
    real(DP), intent(in) :: hkn !< horizontal hydraulic conductivity for cell n (in the direction of cell m)
    real(DP), intent(in) :: hkm !< horizontal hydraulic conductivity for cell m (in the direction of cell n)
    real(DP), intent(in) :: topn !< top of cell n
    real(DP), intent(in) :: topm !< top of cell m
    real(DP), intent(in) :: botn !< bottom of cell n
    real(DP), intent(in) :: botm !< bottom of cell m
    real(DP), intent(in) :: cln !< distance from the center of cell n to the shared face with cell m
    real(DP), intent(in) :: clm !< distance from the center of cell m to the shared face with cell n
    real(DP), intent(in) :: fawidth !< width of cell perpendicular to flow
    !  local
    real(DP) :: thksatn
    real(DP) :: thksatm

    if (ihc == C3D_STAGGERED) then
      thksatn = staggered_thkfrac(topn, botn, satn, topm, botm)
      thksatm = staggered_thkfrac(topm, botm, satm, topn, botn)
    else
      thksatn = satn * (topn - botn)
      thksatm = satm * (topm - botm)
    end if
    condnm = condmean(hkn, hkm, thksatn, thksatm, cln, clm, &
                      fawidth, icellavg)
  end function convertible_standard

  !> @brief Vertical conductance between two cells
  !<
  function vcond(ibdn, ibdm, ictn, ictm, inewton, ivarcv, idewatcv, &
                 condsat, hn, hm, vkn, vkm, satn, satm, topn, topm, botn, &
                 botm, flowarea) result(condnm)
    ! return variable
    real(DP) :: condnm
    ! dummy
    integer(I4B), intent(in) :: ibdn !< cell n active flag
    integer(I4B), intent(in) :: ibdm !< cell m active flag
    integer(I4B), intent(in) :: ictn !< cell n convertible cell flag
    integer(I4B), intent(in) :: ictm !< cell m convertible cell flag
    integer(I4B), intent(in) :: inewton !< flag for Newton-Raphson formulation
    integer(I4B), intent(in) :: ivarcv !< variable vertical conductance flag
    integer(I4B), intent(in) :: idewatcv !< dewatered vertical conductance flag
    real(DP), intent(in) :: condsat !< saturated conductance
    real(DP), intent(in) :: hn !< cell n head
    real(DP), intent(in) :: hm !< cell m head
    real(DP), intent(in) :: vkn !< vertical hydraulic conductivity for cell n (in the direction of cell m)
    real(DP), intent(in) :: vkm !< vertical hydraulic conductivity for cell m (in the direction of cell n)
    real(DP), intent(in) :: satn !< cell n wetted cell fraction
    real(DP), intent(in) :: satm !< cell m wetted cell fraction
    real(DP), intent(in) :: topn !< top of cell n
    real(DP), intent(in) :: topm !< top of cell m
    real(DP), intent(in) :: botn !< bottom of cell n
    real(DP), intent(in) :: botm !< bottom of cell m
    real(DP), intent(in) :: flowarea !< flow area between cell n and m
    ! local
    real(DP) :: satntmp
    real(DP) :: satmtmp
    real(DP) :: bovk1
    real(DP) :: bovk2
    real(DP) :: denom
    !
    ! Either n or m is inactive
    if (ibdn == 0 .or. ibdm == 0) then
      condnm = DZERO
      !
      ! constantcv
    elseif (ivarcv == 0) then
      condnm = condsat
      !
      ! both cells are non-convertible
    elseif (ictn == 0 .and. ictm == 0) then
      condnm = condsat
      !
      ! both cells are fully saturated
    elseif (hn >= topn .and. hm >= topm) then
      condnm = condsat
      !
      ! At least one cell is partially saturated, so recalculate vertical
      ! conductance for this connection
      ! todo: upstream weighting?
    else
      !
      ! Default is for CV correction (dewatered option); use underlying
      !    saturation of 1.
      satntmp = satn
      satmtmp = satm
      if (idewatcv == 0) then
        if (botn > botm) then
          satmtmp = DONE
        else
          satntmp = DONE
        end if
      else
        if (botn > botm) then
          if (satmtmp < DONE) then
            satmtmp = DZERO
          end if
        else
          if (satntmp < DONE) then
            satntmp = DZERO
          end if
        end if
      end if
      bovk1 = satntmp * (topn - botn) * DHALF / vkn
      bovk2 = satmtmp * (topm - botm) * DHALF / vkm
      denom = (bovk1 + bovk2)
      if (denom /= DZERO) then
        condnm = flowarea / denom
      else
        condnm = DZERO
      end if
    end if
  end function vcond

  !> @brief Calculate the conductance between two cells
  !<
  function condmean(k1, k2, thick1, thick2, cl1, cl2, width, iavgmeth)
    ! return variable
    real(DP) :: condmean !< mean conductance between two cells
    ! dummy
    real(DP), intent(in) :: k1 !< hydraulic conductivity for cell n (in the direction of cell m)
    real(DP), intent(in) :: k2 !< hydraulic conductivity for cell m (in the direction of celln)
    real(DP), intent(in) :: thick1 !< saturated thickness for cell 1
    real(DP), intent(in) :: thick2 !< saturated thickness for cell 2
    real(DP), intent(in) :: cl1 !< distance from the center of cell n to the shared face with cell m
    real(DP), intent(in) :: cl2 !< distance from the center of cell m to the shared face with cell n
    real(DP), intent(in) :: width !< width of cell perpendicular to flow
    integer(I4B), intent(in) :: iavgmeth !< averaging method
    ! local
    real(DP) :: t1
    real(DP) :: t2
    real(DP) :: tmean
    real(DP) :: kmean
    real(DP) :: denom
    !
    ! Initialize
    t1 = k1 * thick1
    t2 = k2 * thick2

    ! Averaging method
    select case (iavgmeth)

    case (CCOND_HMEAN)
      if (t1 * t2 > DZERO) then
        condmean = width * t1 * t2 / (t1 * cl2 + t2 * cl1)
      else
        condmean = DZERO
      end if

    case (CCOND_LMEAN)
      if (t1 * t2 > DZERO) then
        tmean = logmean(t1, t2)
      else
        tmean = DZERO
      end if
      condmean = tmean * width / (cl1 + cl2)

    case (CCOND_AMTLMK)
      if (k1 * k2 > DZERO) then
        kmean = logmean(k1, k2)
      else
        kmean = DZERO
      end if
      condmean = kmean * DHALF * (thick1 + thick2) * width / (cl1 + cl2)

    case (CCOND_AMTHMK)
      denom = (k1 * cl2 + k2 * cl1)
      if (denom > DZERO) then
        kmean = k1 * k2 / denom
      else
        kmean = DZERO
      end if
      condmean = kmean * DHALF * (thick1 + thick2) * width
    end select
  end function condmean

  !> @brief Calculate the the logarithmic mean of two double precision numbers
  !!
  !! Use an approximation if the ratio is near 1
  !<
  function logmean(d1, d2)
    ! return variable
    real(DP) :: logmean !< logarithmic mean for two number
    ! dummy
    real(DP), intent(in) :: d1 !< first number
    real(DP), intent(in) :: d2 !< second number
    ! local
    real(DP) :: drat

    drat = d2 / d1
    if (drat <= DLNLOW .or. drat >= DLNHIGH) then
      logmean = (d2 - d1) / log(drat)
    else
      logmean = DHALF * (d1 + d2)
    end if
  end function logmean

  !> @brief Calculate wetted cell thickness at interface between two cells
  !<
  function thksatnm(ibdn, ibdm, ictn, ictm, iupstream, ihc, &
                    hn, hm, satn, satm, topn, topm, botn, botm) result(res)
    ! return variable
    real(DP) :: res !< wetted cell thickness for connection nm
    ! dummy
    integer(I4B), intent(in) :: ibdn !< cell n active flag
    integer(I4B), intent(in) :: ibdm !< cell m active flag
    integer(I4B), intent(in) :: ictn !< cell n convertible cell flag
    integer(I4B), intent(in) :: ictm !< cell m convertible cell flag
    integer(I4B), intent(in) :: iupstream !< flag for upstream weighting
    integer(I4B), intent(in) :: ihc !< connection type
    real(DP), intent(in) :: hn !< cell n head
    real(DP), intent(in) :: hm !< cell m head
    real(DP), intent(in) :: satn !< cell n wetted cell fraction
    real(DP), intent(in) :: satm !< cell m wetted cell fraction
    real(DP), intent(in) :: topn !< top of cell n
    real(DP), intent(in) :: topm !< top of cell m
    real(DP), intent(in) :: botn !< bottom of cell n
    real(DP), intent(in) :: botm !< bottom of cell m
    ! local
    real(DP) :: thksatn
    real(DP) :: thksatm
    real(DP) :: sill_top
    real(DP) :: sill_bot
    !
    ! n or m is inactive
    if (ibdn == 0 .or. ibdm == 0) then
      res = DZERO
      !
      ! both cells are non-convertible
    elseif (ictn == 0 .and. ictm == 0) then
      if (ihc == C3D_STAGGERED) then
        sill_top = min(topn, topm)
        sill_bot = max(botn, botm)

        thksatn = max(sill_top - sill_bot, DZERO)
        thksatm = thksatn
      else
        thksatn = topn - botn
        thksatm = topm - botm
      end if
      res = DHALF * (thksatn + thksatm)
      !
      ! At least one of the cells is convertible
    elseif (iupstream == 1) then
      if (hn > hm) then
        res = satn * (topn - botn)
      else
        res = satm * (topm - botm)
      end if
      !
      ! At least one of the cells is convertible and not upstream weighted
    else
      if (ihc == C3D_STAGGERED) then
        thksatn = staggered_thkfrac(topn, botn, satn, topm, botm)
        thksatm = staggered_thkfrac(topm, botm, satm, topn, botn)
      else
        thksatn = satn * (topn - botn)
        thksatm = satm * (topm - botm)
      end if
      res = DHALF * (thksatn + thksatm)
    end if
  end function thksatnm

  !> @brief Calculate the thickness fraction for staggered grids
  !<
  function staggered_thkfrac(top, bot, sat, topc, botc) result(res)
    ! return variable
    real(DP) :: res !< staggered thickness fraction for cell
    ! dummy
    real(DP) :: top !< top of cell
    real(DP) :: bot !< bottom of cell
    real(DP) :: sat !< cell saturation
    real(DP) :: topc !< top of connected cell
    real(DP) :: botc !< bottom of connected cells
    ! local
    real(DP) :: sill_top
    real(DP) :: sill_bot
    real(DP) :: tp

    sill_top = min(top, topc)
    sill_bot = max(bot, botc)
    tp = bot + sat * (top - bot)
    res = max(min(tp, sill_top) - sill_bot, DZERO)
  end function staggered_thkfrac

end module GwfConductanceUtilsModule
