submodule(SfrModule) SfrModuleSteady
contains

  ! -- documented on the SfrType type-bound procedure in gwf-sfr.f90
  module procedure sfr_calc_steady
  integer(I4B) :: i
  integer(I4B) :: isolve
  integer(I4B) :: iic
  integer(I4B) :: iic2
  integer(I4B) :: iic3
  integer(I4B) :: iic4
  integer(I4B) :: ibflg
  real(DP) :: qmp
  real(DP) :: qsrc
  real(DP) :: qsrcmp
  real(DP) :: tp
  real(DP) :: bt
  real(DP) :: hsfr
  real(DP) :: qc
  real(DP) :: en1
  real(DP) :: en2
  real(DP) :: qen1
  real(DP) :: f1
  real(DP) :: f2
  real(DP) :: qgwf1
  real(DP) :: qgwf2
  real(DP) :: qgwfp
  real(DP) :: qgwfold
  real(DP) :: fhstr1
  real(DP) :: fhstr2
  real(DP) :: d2
  real(DP) :: dpp
  real(DP) :: dx
  real(DP) :: q1
  real(DP) :: q2
  real(DP) :: derv
  real(DP) :: dlh
  real(DP) :: dlhold
  real(DP) :: fp
  real(DP) :: err
  real(DP) :: errold

  d2 = DZERO
  q1 = DZERO
  q2 = DZERO
  qsrc = DZERO
  qgwf = DZERO
  qgwfold = DZERO

  ! calculate the flow at end of the reach
  ! excluding groundwater leakage
  qc = qu + qi + qr - qe + qro + qfrommvr

  ! calculate flow at the middle of the reach
  ! excluding groundwater leakage
  qsrcmp = qu + qi + qfrommvr + DHALF * (qr - qe + qro)
  qmp = qsrcmp ! initial estimate flow at the midpoint

  ! calculate stream depth at the midpoint
  call this%sfr_calc_reach_depth(n, qmp, d1)

  ! calculate sources/sinks for reach
  ! excluding groundwater leakage
  call this%sfr_calc_qsource(n, d1, qsrc)

  ! calculate initial reach stage, downstream flow,
  ! and groundwater leakage
  tp = this%strtop(n)
  bt = tp - this%bthick(n)
  hsfr = d1 + tp
  qd = MAX(qsrc, DZERO)
  qgwf = DZERO
  ! set flag to skip iterations
  isolve = 1
  if (hsfr <= tp .and. hgwf <= tp) isolve = 0
  if (hgwf <= tp .and. qc < DEM30) isolve = 0
  if (this%sfr_gwf_conn(n) == 0) isolve = 0

  ! iterate to achieve solution
  calc_solution: if (isolve /= 0) then

    ! estimate initial end points
    en1 = DZERO
    if (d1 > DEM30) then
    if ((tp - hgwf) > DEM30) then
      en2 = DP9 * d1
    else
      en2 = D1P1 * d1 - (tp - hgwf)
    end if
    else if ((tp - hgwf) > DEM30) then
    en2 = DONE
    else
    en2 = DP99 * (hgwf - tp)
    end if

    ! estimate flow at end points
    ! end point 1
    if (hgwf > tp) then
      call this%sfr_calc_qgwf(n, DZERO, hgwf, qgwf1)
      qgwf1 = -qgwf1
      qen1 = qmp - DHALF * qgwf1
    else
      qgwf1 = DZERO
      qen1 = qsrcmp
    end if
    if (hgwf > bt) then
      call this%sfr_calc_qgwf(n, en2, hgwf, qgwf2)
      qgwf2 = -qgwf2
    else
      call this%sfr_calc_qgwf(n, en2, bt, qgwf2)
      qgwf2 = -qgwf2
    end if
    if (qgwf2 > qsrc) qgwf2 = qsrc

    ! calculate two depths
    call this%sfr_calc_reach_depth(n, (qsrcmp - DHALF * qgwf1), d1)
    call this%sfr_calc_reach_depth(n, (qsrcmp - DHALF * qgwf2), d2)

    ! determine roots
    if (d1 > DEM30) then
      f1 = en1 - d1
    else
      en1 = DZERO
      f1 = en1 - DZERO
    end if
    if (d2 > DEM30) then
      f2 = en2 - d2
      if (f2 < DEM30) en2 = d2
    else
      d2 = DZERO
      f2 = en2 - DZERO
    end if

    ! iterate to find a solution
    dpp = DHALF * (en1 + en2)
    dx = dpp
    iic = 0
    iic2 = 0
    iic3 = 0
    fhstr1 = DZERO
    fhstr2 = DZERO
    qgwfp = DZERO
    dlhold = DZERO
    do i = 1, this%maxsfrit
      ibflg = 0
      d1 = dpp
      d2 = d1 + DTWO * this%deps
      ! calculate q at midpoint at both end points
      call this%sfr_calc_qman(n, d1, q1)
      call this%sfr_calc_qman(n, d2, q2)
      ! calculate groundwater leakage at both end points
      call this%sfr_calc_qgwf(n, d1, hgwf, qgwf1)
      qgwf1 = -qgwf1
      call this%sfr_calc_qgwf(n, d2, hgwf, qgwf2)
      qgwf2 = -qgwf2
      !
      if (qgwf1 >= qsrc) then
        en2 = dpp
        dpp = DHALF * (en1 + en2)
        call this%sfr_calc_qgwf(n, dpp, hgwf, qgwfp)
        qgwfp = -qgwfp
        if (qgwfp > qsrc) qgwfp = qsrc
        call this%sfr_calc_reach_depth(n, (qsrcmp - DHALF * qgwfp), dx)
        ibflg = 1
      else
        fhstr1 = (qsrcmp - DHALF * qgwf1) - q1
        fhstr2 = (qsrcmp - DHALF * qgwf2) - q2
      end if
      !
      if (ibflg == 0) then
        derv = DZERO
        if (abs(d1 - d2) > DZERO) then
          derv = (fhstr1 - fhstr2) / (d1 - d2)
        end if
        if (abs(derv) > DEM30) then
          dlh = -fhstr1 / derv
        else
          dlh = DZERO
        end if
        dpp = d1 + dlh

        ! updated depth outside of endpoints - use bisection instead
        if ((dpp >= en2) .or. (dpp <= en1)) then
        if (abs(dlh) > abs(dlhold) .or. dpp < DEM30) then
          ibflg = 1
          dpp = DHALF * (en1 + en2)
        end if
        end if

        ! check for slow convergence
        ! set flags to determine if the Newton-Raphson method oscillates
        ! or if convergence is slow
        if (qgwf1 * qgwfold < DEM30) then
          iic2 = iic2 + 1
        else
          iic2 = 0
        end if
        if (qgwf1 < DEM30) then
          iic3 = iic3 + 1
        else
          iic3 = 0
        end if
        if (dlh * dlhold < DEM30 .or. ABS(dlh) > ABS(dlhold)) then
          iic = iic + 1
        end if
        iic4 = 0
        if (iic3 > 7 .and. iic > 12) then
          iic4 = 1
        end if

        ! switch to bisection when the Newton-Raphson method oscillates
        ! or when convergence is slow
        if (iic2 > 7 .or. iic > 12 .or. iic4 == 1) then
          ibflg = 1
          dpp = DHALF * (en1 + en2)
        end if
        !
        ! Calculate perturbed gwf flow
        call this%sfr_calc_qgwf(n, dpp, hgwf, qgwfp)
        qgwfp = -qgwfp
        if (qgwfp > qsrc) then
          qgwfp = qsrc
          if (abs(en1 - en2) < this%dmaxchg * DEM6) then
            call this%sfr_calc_reach_depth(n, (qsrcmp - DHALF * qgwfp), dpp)
          end if
        end if
        call this%sfr_calc_reach_depth(n, (qsrcmp - DHALF * qgwfp), dx)
      end if

      ! bisection to update end points
      fp = dpp - dx
      if (ibflg == 1) then
        dlh = fp
        ! change end points
        ! root is between f1 and fp
        if (f1 * fp < DZERO) then
          en2 = dpp
          f2 = fp
          ! root is between fp and f2
        else
          en1 = dpp
          f1 = fp
        end if
        err = min(abs(fp), abs(en2 - en1))
      else
        err = abs(dlh)
      end if

      ! check for convergence and exit if converged
      if (err < this%dmaxchg) then
        d1 = dpp
        qgwf = qgwfp
        qd = qsrc - qgwf
        exit
      end if

      ! save iterates
      errold = err
      dlhold = dlh
      if (ibflg == 1) then
        qgwfold = qgwfp
      else
        qgwfold = qgwf1
      end if

    end do
    ! depth = 0 and hgwf < bt
  else
    call this%sfr_calc_qgwf(n, d1, hgwf, qgwf)
    qgwf = -qgwf

    ! leakage exceeds inflow
    if (qgwf > qsrc) then
      d1 = DZERO
      call this%sfr_calc_qsource(n, d1, qsrc)
      qgwf = qsrc
    end if
    qd = qsrc - qgwf
  end if calc_solution

  end procedure sfr_calc_steady
end submodule
