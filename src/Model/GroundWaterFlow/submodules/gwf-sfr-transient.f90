submodule(SfrModule) SfrModuleTransient
contains

  !> @brief Method for solving for transient reaches
    !!
    !! Method to solve the continuity equation for a SFR package
    !! reach using the kinematic wave approximation.
    !!
  !<
  module procedure sfr_calc_transient
  use TdisModule, only: delt

  integer(I4B) :: igwfconn
  integer(I4B) :: number_picard
  integer(I4B) :: i
  integer(I4B) :: j
  real(DP) :: kinematic_residual
  real(DP) :: kinematic_storage
  real(DP) :: weight
  real(DP) :: celerity
  real(DP) :: courant
  real(DP) :: dq
  real(DP) :: qtol
  real(DP) :: qsrc
  real(DP) :: qlat
  real(DP) :: da
  real(DP) :: db
  real(DP) :: dc
  real(DP) :: dd
  real(DP) :: qa
  real(DP) :: qb
  real(DP) :: qc
  real(DP) :: xsa_a
  real(DP) :: xsa_b
  real(DP) :: xsa_c
  real(DP) :: q
  real(DP) :: q2
  real(DP) :: d
  real(DP) :: d2
  real(DP) :: a
  real(DP) :: a2
  real(DP) :: d1old
  real(DP) :: qd2
  real(DP) :: ad
  real(DP) :: ad2
  real(DP) :: residual
  real(DP) :: residual2
  real(DP) :: residual_final
  real(DP) :: qderv
  real(DP) :: delq
  real(DP) :: delh
  real(DP) :: dd2

  weight = this%storage_weight
  dq = this%deps
  qtol = dq * DTWO

  celerity = DZERO
  qgwf = DZERO

  qlat = qr + qro - qe

  this%usinflow(n) = qu + qi + qfrommvr

  qa = this%usinflowold(n)
  qb = this%dsflowold(n)
  call this%sfr_calc_reach_depth(n, qa, da)
  call this%sfr_calc_reach_depth(n, qb, db)

  qc = this%usinflow(n)
  call this%sfr_calc_reach_depth(n, qc, dc)

  xsa_a = this%calc_area_wet(n, da)
  xsa_b = this%calc_area_wet(n, db)
  xsa_c = this%calc_area_wet(n, dc)

  ! estimate qd
  qd = this%dsflow(n)
  if (qd == DZERO) then
    qd = (qc + qb) * DHALF
  end if
  call this%sfr_calc_reach_depth(n, qd, dd)
  ad = this%calc_area_wet(n, dd)

  ! estimate the depth at the midpoint
  d1 = (dc + dd) * DHALF
  d1old = d1

  ! estimate qgwf
  igwfconn = this%sfr_gwf_conn(n)
  if (igwfconn == 1) then
    q = qu + qi + qr - qe + qro + qfrommvr
    call this%sfr_calc_qgwf(n, d1, hgwf, qgwf)
    qgwf = -qgwf
    if (qgwf > q) then
      qgwf = q
    end if
  end if

  ! calculate maximum wave speed and courant number
  q = qc + qlat - qgwf
  call this%sfr_calc_reach_depth(n, q, d)
  a = this%calc_area_wet(n, d)
  if (d > DZERO) then
    q2 = q + dq
    call this%sfr_calc_reach_depth(n, q2, d2)
    a2 = this%calc_area_wet(n, d2)
    celerity = (q2 - q) / (a2 - a)
    courant = celerity * delt / this%length(n)
    ! write(*,*) this%length(n), courant * delt
  end if

  qlat = qlat / this%length(n)

  number_picard = this%maxsfrpicard
  if (igwfconn == 1) then
    number_picard = this%maxsfrpicard
  else
    number_picard = 1
  end if

  kinematicpicard: do i = 1, number_picard
    if (igwfconn == 1) then
      q = qu + qi + qr - qe + qro + qfrommvr
      call this%sfr_calc_qgwf(n, d1, hgwf, qgwf)
      qgwf = -qgwf
      if (qgwf > q) then
        qgwf = q
      end if
    end if

    qsrc = qlat - qgwf / this%length(n)

    newton: do j = 1, this%maxsfrit
      qd2 = qd + dq
      call this%sfr_calc_reach_depth(n, qd2, dd2)
      ad2 = this%calc_area_wet(n, dd2)

      residual = kinematic_residual(qa, qb, qc, qd, &
                                    xsa_a, xsa_b, xsa_c, ad, &
                                    qsrc, this%length(n), weight, delt, &
                                    courant)

      residual2 = kinematic_residual(qa, qb, qc, qd2, &
                                     xsa_a, xsa_b, xsa_c, ad2, &
                                     qsrc, this%length(n), weight, delt, &
                                     courant)
      qderv = (residual2 - residual) / dq
      if (qderv > DZERO) then
        delq = -residual / qderv
      else
        delq = DZERO
      end if

      if (qd + delq < DEM30) then
        delq = -qd
      end if

      qd = qd + delq

      call this%sfr_calc_reach_depth(n, qd, dd)
      ad = this%calc_area_wet(n, dd)
      residual_final = kinematic_residual(qa, qb, qc, qd, &
                                          xsa_a, xsa_b, xsa_c, ad, &
                                          qsrc, this%length(n), weight, delt, &
                                          courant)

      if (abs(delq) < qtol .and. abs(residual_final) < qtol) then
        exit newton
      end if

    end do newton

    qd = max(qd, DZERO)
    d1 = (dc + dd) * DHALF
    delh = (d1 - d1old)

    if (i > 1 .and. abs(delh) < this%dmaxchg) then
      exit kinematicpicard
    end if

  end do kinematicpicard

  this%storage(n) = kinematic_storage(xsa_a, xsa_b, xsa_c, ad, &
                                      this%length(n), delt, &
                                      courant)

  end procedure sfr_calc_transient

end submodule
