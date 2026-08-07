!> @brief Kinematic wave routing through the unsaturated zone
!<
submodule(UzfCellGroupModule) UzfCellGroupWaves
contains

  !> @brief Reset waves to default values at start of simulation
  !<
  module procedure setwaves
  ! -- local
  real(DP) :: bottom, top
  integer(I4B) :: j
  real(DP) :: thick
  !
  this%flux_to_wt(icell) = DZERO
  this%nwaves(icell) = 1
  thick = this%celtop(icell) - this%water_table(icell)
  do j = 1, this%nwaves_max
    this%wave_depth(j, icell) = DZERO
    this%wave_theta(j, icell) = this%theta_res(icell)
  end do
  !
  ! -- initialize waves for first stress period
  if (thick > DZERO) then
    this%wave_depth(1, icell) = thick
    this%wave_theta(1, icell) = this%theta_init(icell)
    top = this%wave_theta(1, icell) - this%theta_res(icell)
    if (top < DZERO) top = DZERO
    bottom = this%theta_sat(icell) - this%theta_res(icell)
    if (bottom < DZERO) bottom = DZERO
    this%wave_flux(1, icell) = this%vks(icell) &
                               * (top / bottom)**this%bc_eps(icell)
    if (this%wave_theta(1, icell) < this%theta_res(icell)) &
      this%wave_theta(1, icell) = this%theta_res(icell)
    !
    ! -- calculate water stored in the unsaturated zone
    if (top > DZERO) then
      this%wave_speed(1, icell) = DZERO
    else
      this%wave_flux(1, icell) = DZERO
      this%wave_speed(1, icell) = DZERO
    end if
    !
    ! no unsaturated zone
  else
    this%wave_flux(1, icell) = DZERO
    this%wave_depth(1, icell) = DZERO
    this%wave_speed(1, icell) = DZERO
    this%wave_theta(1, icell) = this%theta_res(icell)
  end if
  end procedure setwaves

  !> @brief Prepare and route waves over time step
  !<
  module procedure routewaves
  ! -- local
  real(DP) :: thick, thickold
  integer(I4B) :: idelt, j, ik
  !
  this%flux_to_wt(icell) = DZERO
  this%et_uz(icell) = DZERO
  thick = this%celtop(icell) - this%water_table(icell)
  thickold = this%celtop(icell) - this%water_table_old(icell)
  !
  ! -- no uz, clear waves
  if (thickold < DZERO) then
    do j = 1, 6
      this%wave_theta(j, icell) = this%theta_res(icell)
      this%wave_depth(j, icell) = DZERO
      this%wave_speed(j, icell) = DZERO
      this%wave_flux(j, icell) = DZERO
      this%nwaves(icell) = 1
    end do
  end if
  idelt = 1
  do ik = 1, idelt
    call this%uzflow(thick, thickold, delt, ietflag, icell, ierr)
    if (ierr > 0) return
    totfluxtot = totfluxtot + this%flux_to_wt(icell)
  end do
  end procedure routewaves

  !> @brief Move waves within a cell by shft positions
  !!
  !! Waves strt to stp are overwritten by the waves shft positions away.
  !! cntr is the loop increment; it must run the loop in the direction that
  !! keeps the source ahead of the destination.
  !<
  module procedure shift_waves
  ! -- local
  integer(I4B) :: j
  !
  do j = strt, stp, cntr
    this%wave_theta(j, icell) = this%wave_theta(j + shft, icell)
    this%wave_depth(j, icell) = this%wave_depth(j + shft, icell)
    this%wave_flux(j, icell) = this%wave_flux(j + shft, icell)
    this%wave_speed(j, icell) = this%wave_speed(j + shft, icell)
  end do
  end procedure shift_waves

  !> @brief Save the wave train of one cell into a snapshot
  !<
  module procedure store_waves
  ! -- local
  integer(I4B) :: j
  !
  do j = 1, this%nwaves(icell)
    store%theta(j) = this%wave_theta(j, icell)
    store%depth(j) = this%wave_depth(j, icell)
    store%flux(j) = this%wave_flux(j, icell)
    store%speed(j) = this%wave_speed(j, icell)
  end do
  store%nwaves = this%nwaves(icell)
  end procedure store_waves

  !> @brief Restore the wave train of one cell from a snapshot
  !<
  module procedure load_waves
  ! -- local
  integer(I4B) :: j
  !
  do j = 1, store%nwaves
    this%wave_theta(j, icell) = store%theta(j)
    this%wave_depth(j, icell) = store%depth(j)
    this%wave_flux(j, icell) = store%flux(j)
    this%wave_speed(j, icell) = store%speed(j)
  end do
  this%nwaves(icell) = store%nwaves
  end procedure load_waves

  !> @brief Method of Characteristics solution for kinematic wave equation
  !<
  module procedure uzflow
  ! -- local
  real(DP) :: dflux_surf, time, eps_thick, eps_flux
  real(DP) :: thetadif, thetab, fluxb, oldsflx
  integer(I4B) :: trail_added, single_wave
  !
  time = DZERO
  this%flux_to_wt(icell) = DZERO
  trail_added = 0
  oldsflx = this%wave_flux(this%nwaves(icell), icell)
  call factors(eps_thick, eps_flux)
  !
  ! -- check for falling or rising water table
  if ((thick - thickold) > eps_thick) then
    thetadif = abs(this%wave_theta(1, icell) - this%theta_res(icell))
    if (thetadif > DEM6) then
      call this%shift_waves(icell, -1, &
                            this%nwaves(icell) + 1, 2, -1)
      if (this%wave_depth(2, icell) < DEM30) &
        this%wave_depth(2, icell) = (this%ntrailwaves + DTWO) * DEM6
      if (this%wave_theta(2, icell) > this%theta_res(icell)) then
        this%wave_speed(2, icell) = &
          this%wave_flux(2, icell) / &
          (this%wave_theta(2, icell) - this%theta_res(icell))
      else
        this%wave_speed(2, icell) = DZERO
      end if
      this%wave_theta(1, icell) = this%theta_res(icell)
      this%wave_flux(1, icell) = DZERO
      this%wave_speed(1, icell) = DZERO
      this%wave_depth(1, icell) = thick
      this%nwaves(icell) = this%nwaves(icell) + 1
      if (this%nwaves(icell) >= this%nwaves_max) then
        ! -- too many waves error
        ierr = 1
        return
      end if
    else
      this%wave_depth(1, icell) = thick
    end if
  end if
  thetab = this%wave_theta(1, icell)
  fluxb = this%wave_flux(1, icell)
  this%flux_to_wt(icell) = DZERO
  single_wave = 0
  dflux_surf = (this%surf_infil(icell) - this%wave_flux(this%nwaves(icell), &
                                                        icell))
  !
  ! -- increase new waves in infiltration changes
  if (dflux_surf > eps_flux .OR. dflux_surf < -eps_flux) then
    this%nwaves(icell) = this%nwaves(icell) + 1
    if (this%nwaves(icell) >= this%nwaves_max) then
      !
      ! -- too many waves error
      ierr = 1
      return
    end if
  else if (this%nwaves(icell) == 1) then
    single_wave = 1
  end if
  if (this%nwaves(icell) > 1) then
    if (dflux_surf < -eps_flux) then
      call this%trailwav(icell, ierr)
      if (ierr > 0) return
      trail_added = 1
    end if
    call this%leadwav(time, single_wave, trail_added, thetab, fluxb, dflux_surf, &
                      eps_flux, delt, icell)
  end if
  if (single_wave == 1) then
    this%flux_to_wt(icell) = this%flux_to_wt(icell) + &
                             (delt - time) * this%wave_flux(1, icell)
    time = DZERO
    single_wave = 0
  end if
  !
  ! -- simulate et
  if (ietflag > 0) call this%uzet(icell, delt, ietflag, ierr)
  if (ierr > 0) return
  end procedure uzflow

  !> @brief Calculate unit specific tolerances
  !<
  module procedure factors
  ! -- local
  real(DP) :: factor1
  real(DP) :: factor2
  !
  factor1 = DONE
  factor2 = DONE
  eps_thick = DEM9
  eps_flux = DEM9
  if (ITMUNI == 1) then
    factor1 = DONE / 86400.D0
  else if (ITMUNI == 2) then
    factor1 = DONE / 1440.D0
  else if (ITMUNI == 3) then
    factor1 = DONE / 24.0D0
  else if (ITMUNI == 5) then
    factor1 = 365.0D0
  end if
  factor2 = DONE / 0.3048
  eps_thick = eps_thick * factor1 * factor2
  eps_flux = eps_flux * factor1 * factor2
  end procedure factors

  !> @brief Create and set trail waves
  !<
  module procedure trailwav
  ! -- local
  real(DP) :: theta_surf, theta_step, ftrail, eps_m1
  real(DP) :: dtheta_inv
  real(DP) :: sat
  real(DP) :: flux1, flux2, theta1, theta2
  real(DP) :: fnuminc
  integer(I4B) :: j, jj, jk, nwaves_m1
  !
  eps_m1 = dble(this%bc_eps(icell)) - DONE
  dtheta_inv = DONE / (this%theta_sat(icell) - this%theta_res(icell))
  nwaves_m1 = this%nwaves(icell) - 1
  !
  ! -- initialize trailwaves
  sat = (this%surf_infil(icell) / this%vks(icell))**(DONE / this%bc_eps(icell))
  theta_surf = (sat * (this%theta_sat(icell) - this%theta_res(icell))) + &
               this%theta_res(icell)
  if (this%wave_theta(nwaves_m1, icell) - theta_surf > DEM9) then
    fnuminc = DZERO
    do jk = 1, this%ntrailwaves
      fnuminc = fnuminc + float(jk)
    end do
    theta_step = (this%wave_theta(nwaves_m1, &
                                  icell) - theta_surf) / (fnuminc - DONE)
    jj = this%ntrailwaves
    ftrail = dble(this%ntrailwaves) + DONE
    do j = this%nwaves(icell), this%nwaves(icell) + this%ntrailwaves - 1
      if (j > this%nwaves_max) then
        ! -- too many waves error
        ierr = 1
        return
      end if
      if (j > this%nwaves(icell)) then
        this%wave_theta(j, icell) = this%wave_theta(j - 1, icell) &
                                    - ((ftrail - float(jj)) * theta_step)
      else
        this%wave_theta(j, icell) = this%wave_theta(j - 1, icell) - DEM9
      end if
      jj = jj - 1
      if (this%wave_theta(j, icell) <= this%theta_res(icell) + DEM9) &
        this%wave_theta(j, icell) = this%theta_res(icell) + DEM9
      sat = (this%wave_theta(j, icell) - this%theta_res(icell)) * dtheta_inv
      this%wave_flux(j, icell) = this%vks(icell) * (sat**this%bc_eps(icell))
      theta2 = this%wave_theta(j - 1, icell)
      flux2 = this%wave_flux(j - 1, icell)
      flux1 = this%wave_flux(j, icell)
      theta1 = this%wave_theta(j, icell)
      this%wave_speed(j, icell) = &
        leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                  this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
      this%wave_depth(j, icell) = DZERO
      if (j == this%nwaves(icell)) then
        this%wave_depth(j, icell) = this%wave_depth(j, icell) + &
                                    (this%ntrailwaves + 1) * DEM9
      else
        this%wave_depth(j, icell) = this%wave_depth(j - 1, icell) - DEM9
      end if
    end do
    this%nwaves(icell) = this%nwaves(icell) + this%ntrailwaves - 1
    if (this%nwaves(icell) >= this%nwaves_max) then
      ! -- too many waves error
      ierr = 1
      return
    end if
  else
    this%wave_depth(this%nwaves, icell) = DZERO
    this%wave_flux(this%nwaves, icell) = &
      this%vks(icell) * (((this%wave_theta(this%nwaves, icell) - &
                           this%theta_res(icell)) * dtheta_inv)** &
                         this%bc_eps(icell))
    this%wave_theta(this%nwaves, icell) = theta_surf
    theta2 = this%wave_theta(this%nwaves(icell) - 1, icell)
    flux2 = this%wave_flux(this%nwaves(icell) - 1, icell)
    flux1 = this%wave_flux(this%nwaves(icell), icell)
    theta1 = this%wave_theta(this%nwaves(icell), icell)
    this%wave_speed(this%nwaves(icell), icell) = &
      leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
  end if
  end procedure trailwav

  !> @brief Create a lead wave and route over time step
  !<
  module procedure leadwav
  ! -- local
  real(DP) :: dt_to_bottom, dt_first, fcheck
  real(DP) :: eps_m1, time_new, bottom, time_left
  real(DP) :: dtheta_inv, diff, flux_bot_prev
  real(DP) :: flux1, flux2, theta1, theta2, ftest
  integer(I4B) :: wave_exited, iremove, j, jmerge
  integer(I4B) :: nwavp1, j_first
  !
  ftest = DZERO
  eps_m1 = dble(this%bc_eps(icell)) - DONE
  dtheta_inv = DONE / (this%theta_sat(icell) - this%theta_res(icell))
  !
  ! -- initialize new wave
  if (trail_added == 0) then
    if (dflux_surf > eps_flux) then
      this%wave_flux(this%nwaves(icell), icell) = this%surf_infil(icell)
      if (this%wave_flux(this%nwaves(icell), icell) < DEM30) &
        this%wave_flux(this%nwaves(icell), icell) = DZERO
      this%wave_theta(this%nwaves(icell), icell) = &
        (((this%wave_flux(this%nwaves(icell), icell) / this%vks(icell))** &
          (DONE / this%bc_eps(icell))) * &
         (this%theta_sat(icell) - this%theta_res(icell))) + &
        this%theta_res(icell)
      theta2 = this%wave_theta(this%nwaves(icell), icell)
      flux2 = this%wave_flux(this%nwaves(icell), icell)
      flux1 = this%wave_flux(this%nwaves(icell) - 1, icell)
      theta1 = this%wave_theta(this%nwaves(icell) - 1, icell)
      this%wave_speed(this%nwaves(icell), icell) = &
        leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                  this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
      this%wave_depth(this%nwaves(icell), icell) = DZERO
    end if
  end if
  !
  ! -- route all waves and interception of waves over times step
  diff = DONE
  time_left = DZERO
  wave_exited = 0
  flux_bot_prev = this%wave_flux(1, icell)
  if (this%nwaves(icell) == 0) single_wave = 1
  if (single_wave /= 1) then
    do while (diff > DEM6)
      nwavp1 = this%nwaves(icell) + 1
      time_left = delt - Time
      do j = 1, this%nwaves(icell)
        this%overtake_time(j) = DEP20
        this%wave_merges(j) = 0
      end do
      dt_first = time_left
      if (this%nwaves(icell) > 2) then
        j = 2
        !
        ! -- calculate time until wave overtakes wave ahead
        nwavp1 = this%nwaves(icell) + 1
        do while (j < nwavp1)
          ftest = this%wave_speed(j - 1, icell) - this%wave_speed(j, icell)
          if (abs(ftest) > DEM30) then
            this%overtake_time(j) = (this%wave_depth(j, icell) - &
                                     this%wave_depth(j - 1, icell)) / ftest
            if (this%overtake_time(j) < DEM30) this%overtake_time(j) = DEP20
          end if
          j = j + 1
        end do
      end if
      !
      ! - calc time until wave reaches bottom of cell
      dt_to_bottom = DEP20
      if (this%nwaves(icell) > 1) then
        if (this%wave_speed(2, icell) > DZERO) then
          bottom = this%wave_speed(2, icell)
          if (bottom < DEM15) bottom = DEM15
          dt_to_bottom = (this%wave_depth(1, icell) - &
                          this%wave_depth(2, icell)) / bottom
          if (dt_to_bottom < DZERO) dt_to_bottom = DEM12
        end if
      end if
      !
      ! -- calc time for wave interception
      j_first = 0
      do j = this%nwaves(icell), 3, -1
        if (dt_first - this%overtake_time(j) > -DEM9) then
          this%wave_merges(j) = 1
          j_first = j
          dt_first = this%overtake_time(j)
        end if
      end do
      do j = 3, this%nwaves(icell)
        if (dt_first - this%overtake_time(j) < DEM9) then
          if (j /= j_first) this%wave_merges(j) = 0
        end if
      end do
      !
      ! -- what happens first, waves hits bottom or interception
      iremove = 0
      time_new = Time
      fcheck = (Time + dt_first) - delt
      if (dt_first < DEM7) fcheck = -DONE
      if (dt_to_bottom < dt_first .AND. Time + dt_to_bottom < delt) then
        j = 2
        do while (j < nwavp1)
          !
          ! -- route waves
          this%wave_depth(j, icell) = this%wave_depth(j, icell) + &
                                      this%wave_speed(j, icell) * dt_to_bottom
          j = j + 1
        end do
        fluxb = this%wave_flux(2, icell)
        thetab = this%wave_theta(2, icell)
        wave_exited = 1
        call this%shift_waves(icell, 1, 1, &
                              this%nwaves(icell) - 1, 1)
        iremove = 1
        time_new = time + dt_to_bottom
        this%wave_speed(1, icell) = DZERO
        !
        ! -- do waves intercept before end of time step
      else if (fcheck < DZERO .AND. this%nwaves(icell) > 2) then
        j = 2
        do while (j < nwavp1)
          this%wave_depth(j, icell) = this%wave_depth(j, icell) + &
                                      this%wave_speed(j, icell) * dt_first
          j = j + 1
        end do
        !
        ! -- combine waves that intercept, remove a wave
        j = 3
        jmerge = j
        do while (j < this%nwaves(icell) + 1)
          if (this%wave_merges(j) == 1) then
            jmerge = j
            theta2 = this%wave_theta(j, icell)
            flux2 = this%wave_flux(j, icell)
            if (j == 3) then
              flux1 = fluxb
              theta1 = thetab
            else
              flux1 = this%wave_flux(j - 2, icell)
              theta1 = this%wave_theta(j - 2, icell)
            end if
            this%wave_speed(j, icell) = &
              leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                        this%theta_res(icell), this%bc_eps(icell), &
                        this%vks(icell))
            !
            ! -- update waves
            call this%shift_waves(icell, 1, jmerge - 1, &
                                  this%nwaves(icell) - 1, 1)
            jmerge = this%nwaves(icell) + 1
            iremove = iremove + 1
          end if
          j = j + 1
        end do
        time_new = time_new + dt_first
        !
        ! -- calc. total flux to bottom during remaining time in step
      else
        j = 2
        do while (j < nwavp1)
          this%wave_depth(j, icell) = this%wave_depth(j, icell) + &
                                      this%wave_speed(j, icell) * time_left
          j = j + 1
        end do
        time_new = delt
      end if
      this%flux_to_wt(icell) = this%flux_to_wt(icell) + flux_bot_prev &
                               * (time_new - time)
      if (wave_exited == 1) then
        flux_bot_prev = this%wave_flux(1, icell)
        wave_exited = 0
      end if
      !
      ! -- remove dead waves
      this%nwaves(icell) = this%nwaves(icell) - iremove
      time = time_new
      diff = delt - Time
      if (this%nwaves(icell) == 1) then
        single_wave = 1
        exit
      end if
    end do
  end if
  end procedure leadwav

  !> @brief Sums up mobile water over depth interval
  !<
  module procedure unsat_stor
  ! -- local
  real(DP) :: fm
  integer(I4B) :: jabove, j, nwaves_m1
  !
  fm = DZERO
  jabove = this%nwaves(icell) + 1
  nwaves_m1 = this%nwaves(icell) - 1
  if (d1 > this%wave_depth(1, icell)) d1 = this%wave_depth(1, icell)
  !
  ! -- jabove is the deepest wave lying above depth d1. Wave depth decreases
  !    with index, so the first wave that qualifies is the one wanted.
  do j = 1, this%nwaves(icell)
    if (this%wave_depth(j, icell) - d1 < -DEM30) then
      jabove = j
      exit
    end if
  end do
  if (jabove > this%nwaves(icell)) then
    fm = fm + (this%wave_theta(this%nwaves(icell), &
                               icell) - this%theta_res(icell)) * d1
  elseif (this%nwaves(icell) > 1) then
    if (jabove > 1) then
      fm = fm + (this%wave_theta(jabove - 1, icell) - this%theta_res(icell)) &
           * (d1 - this%wave_depth(jabove, icell))
    end if
    do j = jabove, nwaves_m1
      fm = fm + (this%wave_theta(j, icell) - this%theta_res(icell)) &
           * (this%wave_depth(j, icell) &
              - this%wave_depth(j + 1, icell))
    end do
    fm = fm + (this%wave_theta(this%nwaves(icell), icell) - &
               this%theta_res(icell)) * &
         this%wave_depth(this%nwaves(icell), icell)
  else
    fm = fm + (this%wave_theta(1, icell) - this%theta_res(icell)) * d1
  end if
  unsat_stor = fm
  end procedure unsat_stor

  !> @brief Update to new state of uz at end of time step
  !<
  module procedure update_wav
  ! -- local
  real(DP) :: bot, depthsave, top
  real(DP) :: thick, dtheta_inv
  integer(I4B) :: nwaves_hold, j, jabove
  !
  bot = this%water_table(icell)
  top = this%celtop(icell)
  thick = top - bot
  nwaves_hold = this%nwaves(icell)
  if (itest == 1) then
    this%wave_flux(1, icell) = DZERO
    this%wave_theta(1, icell) = this%theta_res(icell)
    return
  end if
  if (iss == 1) then
    if (this%theta_sat(icell) - this%theta_res(icell) < DEM7) then
      dtheta_inv = DONE / DEM7
    else
      dtheta_inv = DONE / (this%theta_sat(icell) - this%theta_res(icell))
    end if
    this%flux_to_wt(icell) = this%surf_infil(icell) * delt
    this%water_table_old(icell) = this%water_table(icell)
    this%wave_theta(1, icell) = this%theta_init(icell)
    this%wave_flux(1, icell) = &
      this%vks(icell) * (((this%wave_theta(1, icell) - this%theta_res(icell)) &
                          * dtheta_inv)**this%bc_eps(icell))
    this%wave_depth(1, icell) = thick
    this%wave_speed(1, icell) = thick
    this%nwaves(icell) = 1
  else
    !
    ! -- water table rises through waves
    if (this%water_table(icell) - this%water_table_old(icell) > DEM30) then
      depthsave = this%wave_depth(1, icell)
      jabove = 0
      do j = 1, this%nwaves(icell)
        if (this%wave_depth(j, icell) - thick < -DEM30) then
          jabove = j
          exit
        end if
      end do
      this%wave_depth(1, icell) = thick
      if (jabove > 1) then
        this%wave_speed(1, icell) = DZERO
        this%nwaves(icell) = this%nwaves(icell) - jabove + 2
        this%wave_theta(1, icell) = this%wave_theta(jabove - 1, icell)
        this%wave_flux(1, icell) = this%wave_flux(jabove - 1, icell)
        if (jabove > 2) call this%shift_waves(icell, jabove - 2, 2, &
                                              nwaves_hold - (jabove - 2), 1)
      elseif (jabove == 0) then
        this%wave_speed(1, icell) = DZERO
        this%wave_theta(1, icell) = this%wave_theta(this%nwaves(icell), icell)
        this%wave_flux(1, icell) = this%wave_flux(this%nwaves(icell), icell)
        this%nwaves(icell) = 1
      end if
    end if
    !
    ! -- calculate new unsat. storage
    if (thick <= DZERO) then
      this%wave_speed(1, icell) = DZERO
      this%nwaves(icell) = 1
      this%wave_theta(1, icell) = this%theta_res(icell)
      this%wave_flux(1, icell) = DZERO
    end if
    this%water_table_old(icell) = this%water_table(icell)
  end if
  end procedure update_wav

  !> @brief Calculate recharge due to a rise in the gwf head
  !<
  module procedure uz_rise
  ! -- local
  real(DP) :: fm1, fm2, d1
  !
  if (this%water_table(icell) - this%water_table_old(icell) > DEM30) then
    d1 = this%celtop(icell) - this%water_table_old(icell)
    fm1 = this%unsat_stor(icell, d1)
    d1 = this%celtop(icell) - this%water_table(icell)
    fm2 = this%unsat_stor(icell, d1)
    totfluxtot = totfluxtot + (fm1 - fm2)
  end if
  end procedure uz_rise

  !> @brief Determine the water content at a specific depth
  !!
  !! Because UZF-calculated waves are internal to UZF objects, different water
  !! contents exists at different depths.
  !<
  module procedure get_water_content_at_depth
  ! -- local
  real(DP) :: d1
  real(DP) :: d2
  real(DP) :: f1
  real(DP) :: f2
  !
  if (this%water_table(icell) < this%celtop(icell)) then
    if (this%celtop(icell) - depth > this%water_table(icell)) then
      d1 = depth - DEM3
      d2 = depth + DEM3
      f1 = this%unsat_stor(icell, d1)
      f2 = this%unsat_stor(icell, d2)
      theta_at_depth = this%theta_res(icell) + (f2 - f1) / (d2 - d1)
    else
      theta_at_depth = this%theta_sat(icell)
    end if
  else
    theta_at_depth = this%theta_sat(icell)
  end if
  end procedure get_water_content_at_depth

  !> @brief Calculate and return the cell-based water content value
  !<
  module procedure get_wcnew
  ! -- local
  real(DP) :: top
  real(DP) :: bot
  real(DP) :: theta_r
  real(DP) :: thk
  real(DP) :: hgwf
  real(DP) :: fm
  real(DP) :: d
  !
  hgwf = this%water_table(icell)
  top = this%celtop(icell)
  bot = this%celbot(icell)
  thk = top - max(bot, hgwf)
  if (thk > DZERO) then
    theta_r = this%theta_res(icell)
    d = thk
    fm = this%unsat_stor(icell, d)
    watercontent = fm / thk
    watercontent = watercontent + theta_r
  else
    watercontent = DZERO
  end if
  end procedure get_wcnew

  !> @brief Calculates waves speed from dflux/dtheta
  !<
  module procedure leadspeed
  ! -- local
  real(DP) :: comp1, comp2, thsrinv, epsfksths
  real(DP) :: eps_m1, fhold, comp3
  !
  eps_m1 = eps - DONE
  thsrinv = DONE / (thts - thtr)
  epsfksths = eps * vks * thsrinv
  comp1 = theta2 - theta1
  comp2 = abs(flux2 - flux1)
  comp3 = theta1 - thtr
  if (comp2 < DEM15) flux2 = flux1 + DEM15
  if (abs(comp1) < DEM30) then
    if (comp3 > DEM30) fhold = (comp3 * thsrinv)**eps
    if (fhold < DEM30) fhold = DEM30
    leadspeed = epsfksths * (fhold**eps_m1)
  else
    leadspeed = (flux2 - flux1) / (theta2 - theta1)
  end if
  if (leadspeed < DEM30) leadspeed = DEM30
  end procedure leadspeed

end submodule UzfCellGroupWaves
