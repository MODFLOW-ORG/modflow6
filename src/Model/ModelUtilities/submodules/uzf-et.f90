!> @brief Evapotranspiration from the unsaturated and saturated zones
!<
submodule(UzfCellGroupModule) UzfCellGroupEt
contains

  !> @brief Subtract aet from pet to calculate residual et for gw
  !<
  module procedure setgwpet
  ! -- modules
  use TdisModule, only: delt
  ! -- local
  real(DP) :: pet
  !
  pet = DZERO
  !
  ! -- reduce pet for gw by uzet
  pet = this%pet(icell) - this%et_uz(icell) / delt
  if (pet < DZERO) pet = DZERO
  this%gw_pet(icell) = pet
  end procedure setgwpet

  !> @brief Subtract aet from pet to calculate residual et for deeper cells
  !<
  module procedure setbelowpet
  ! -- modules
  use TdisModule, only: delt
  ! -- local
  real(DP) :: pet
  !
  pet = DZERO
  !
  ! -- transfer unmet pet to lower cell
  !
  if (this%ext_depth_uz(jbelow) > DEM3) then
    pet = this%pet(icell) - this%et_uz(icell) / delt - &
          this%gwet(icell) / this%uzfarea(icell)
    if (pet < DZERO) pet = DZERO
  end if
  this%pet(jbelow) = pet
  end procedure setbelowpet

  !> @brief Calculate gwf et using residual uzf pet
  !<
  module procedure simgwet
  ! -- local
  real(DP) :: s, x, c, b, et
  !
  this%gwet(icell) = DZERO
  trhs = DZERO
  thcof = DZERO
  det = DZERO
  s = this%landtop(icell)
  x = this%ext_depth(icell)
  c = this%gw_pet(icell)
  b = this%celbot(icell)
  if (b > hgwf) return
  if (x < DEM6) return
  if (igwetflag == 1) then
    et = etfunc_lin(s, x, c, det, trhs, thcof, hgwf, &
                    this%celtop(icell), this%celbot(icell))
  else if (igwetflag == 2) then
    et = etfunc_nlin(s, x, c, det, trhs, thcof, hgwf)
  end if
  ! this%gwet(icell) = et * this%uzfarea(icell)
  trhs = trhs * this%uzfarea(icell)
  thcof = thcof * this%uzfarea(icell)
  this%gwet(icell) = trhs - (thcof * hgwf)
  ! write(99,*)'in group', icell, this%gwet(icell)
  end procedure simgwet

  !> @brief Remove water from uz due to et
  !<
  module procedure uzet
  ! -- local
  real(DP) :: diff
  real(DP) :: thetaout
  real(DP) :: fm
  real(DP) :: st
  real(DP) :: dtheta_inv
  real(DP) :: epsfksthts
  real(DP) :: fmp
  real(DP) :: fktho
  real(DP) :: theta1
  real(DP) :: theta2
  real(DP) :: flux1
  real(DP) :: flux2
  real(DP) :: hcap
  real(DP) :: factor
  real(DP) :: tho
  real(DP) :: depth
  real(DP) :: extwc1
  real(DP) :: theta_min
  real(DP) :: sat
  real(DP) :: petsub
  integer(I4B) :: j
  integer(I4B) :: jext
  integer(I4B) :: jwet
  integer(I4B) :: found
  integer(I4B) :: numadd
  integer(I4B) :: k
  integer(I4B) :: nwv
  integer(I4B) :: itest
  !
  this%et_uz(icell) = DZERO
  if (this%ext_depth_uz(icell) < DEM7) return
  petsub = this%root_act(icell) * this%pet(icell) * &
           this%ext_depth_uz(icell) / this%ext_depth(icell)
  thetaout = delt * petsub / this%ext_depth(icell)
  if (ietflag == 1) thetaout = delt * this%pet(icell) / this%ext_depth(icell)
  if (thetaout < DEM10) return
  depth = this%wave_depth(1, icell)
  st = this%unsat_stor(icell, depth)
  if (st < DEM4) return
  !
  ! -- store original wave characteristics so the aet-to-pet loop can retry
  nwv = this%nwaves(icell)
  itest = 0
  call this%store_waves(icell, this%etsav)
  factor = DONE
  this%et_uz(icell) = DZERO
  if (this%theta_sat(icell) - this%theta_res(icell) < DEM7) then
    dtheta_inv = 1.0 / DEM7
  else
    dtheta_inv = DONE / (this%theta_sat(icell) - this%theta_res(icell))
  end if
  epsfksthts = this%bc_eps(icell) * this%vks(icell) * dtheta_inv
  this%et_uz(icell) = DZERO
  fmp = DZERO
  extwc1 = this%theta_ext(icell) - this%theta_res(icell)
  if (extwc1 < DEM6) extwc1 = DEM7
  !
  ! -- lowest water content et can draw a wave down to. Formed from extwc1
  !    rather than from theta_ext directly, because extwc1 is clamped above.
  theta_min = this%theta_res(icell) + extwc1
  numadd = 0
  fm = st
  k = 0
  !
  ! -- loop for reducing aet to pet when et is head dependent
  do while (itest == 0)
    k = k + 1
    if (k > 1 .AND. ABS(fmp - petsub) > DEM5 * petsub) then
      factor = factor / (fm / petsub)
    end if
    !
    ! -- one wave shallower than extdp
    if (this%nwaves(icell) == 1 .AND. &
        this%wave_depth(1, icell) <= this%ext_depth_uz(icell)) then
      if (ietflag == 2) then
        tho = this%wave_theta(1, icell)
        fktho = this%wave_flux(1, icell)
        hcap = this%caph(icell, tho)
        thetaout = this%rate_et_z(icell, factor, fktho, hcap)
      end if
      if ((this%wave_theta(1, icell) - thetaout) > theta_min) then
        this%wave_theta(1, icell) = this%wave_theta(1, icell) - thetaout
        sat = (this%wave_theta(1, icell) - this%theta_res(icell)) * dtheta_inv
        this%wave_flux(1, icell) = this%vks(icell) * (sat**this%bc_eps(icell))
      else if (this%wave_theta(1, icell) > theta_min) then
        this%wave_theta(1, icell) = theta_min
        sat = (this%wave_theta(1, icell) - this%theta_res(icell)) * dtheta_inv
        this%wave_flux(1, icell) = this%vks(icell) * (sat**this%bc_eps(icell))
      end if
      !
      ! -- all waves shallower than extinction depth
    else if (this%nwaves(icell) > 1 .AND. &
             this%wave_depth(this%nwaves(icell), &
                             icell) > this%ext_depth_uz(icell)) then
      if (ietflag == 2) then
        tho = this%wave_theta(this%nwaves(icell), icell)
        fktho = this%wave_flux(this%nwaves(icell), icell)
        hcap = this%caph(icell, tho)
        thetaout = this%rate_et_z(icell, factor, fktho, hcap)
      end if
      if (this%nwaves(icell) + 1 > this%nwaves_max) then
        !
        ! -- too many waves error
        ierr = 1
        goto 500
      end if
      if (this%wave_theta(this%nwaves(icell), icell) - thetaout > &
          theta_min) then
        this%wave_theta(this%nwaves(icell) + 1, icell) = &
          this%wave_theta(this%nwaves(icell), icell) - thetaout
        numadd = 1
      else if (this%wave_theta(this%nwaves(icell), icell) > &
               theta_min) then
        this%wave_theta(this%nwaves(icell) + 1, icell) = theta_min
        numadd = 1
      end if
      if (numadd == 1) then
        this%wave_flux(this%nwaves(icell) + 1, icell) = &
          this%vks(icell) * &
          (((this%wave_theta(this%nwaves(icell) + 1, icell) - &
             this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
        theta2 = this%wave_theta(this%nwaves(icell) + 1, icell)
        flux2 = this%wave_flux(this%nwaves(icell) + 1, icell)
        flux1 = this%wave_flux(this%nwaves(icell), icell)
        theta1 = this%wave_theta(this%nwaves(icell), icell)
        this%wave_speed(this%nwaves(icell) + 1, icell) = &
          leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                    this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
        this%wave_depth(this%nwaves(icell) + 1, icell) = this%ext_depth_uz(icell)
        this%nwaves(icell) = this%nwaves(icell) + 1
        if (this%nwaves(icell) > this%nwaves_max) then
          !
          ! -- too many waves error
          ierr = 1
          goto 500
        end if
      else
        numadd = 0
      end if
      !
      ! -- one wave below extinction depth
    else if (this%nwaves(icell) == 1) then
      if (this%nwaves(icell) + 1 > this%nwaves_max) then
        !
        ! -- too many waves error
        ierr = 1
        goto 500
      end if
      if (ietflag == 2) then
        tho = this%wave_theta(1, icell)
        fktho = this%wave_flux(1, icell)
        hcap = this%caph(icell, tho)
        thetaout = this%rate_et_z(icell, factor, fktho, hcap)
      end if
      if ((this%wave_theta(1, icell) - thetaout) > theta_min) then
        if (thetaout > DEM30) then
          this%wave_theta(2, icell) = this%wave_theta(1, icell) - thetaout
          sat = (this%wave_theta(2, icell) - this%theta_res(icell)) * dtheta_inv
          this%wave_flux(2, icell) = this%vks(icell) * (sat**this%bc_eps(icell))
          this%wave_depth(2, icell) = this%ext_depth_uz(icell)
          theta2 = this%wave_theta(2, icell)
          flux2 = this%wave_flux(2, icell)
          flux1 = this%wave_flux(1, icell)
          theta1 = this%wave_theta(1, icell)
          this%wave_speed(2, icell) = &
            leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                      this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
          this%nwaves(icell) = this%nwaves(icell) + 1
          if (this%nwaves(icell) > this%nwaves_max) then
            !
            ! -- too many waves error
            ierr = 1
            goto 500
          end if
        end if
      else if (this%wave_theta(1, icell) > theta_min) then
        if (thetaout > DEM30) then
          this%wave_theta(2, icell) = theta_min
          sat = (this%wave_theta(2, icell) - this%theta_res(icell)) * dtheta_inv
          this%wave_flux(2, icell) = this%vks(icell) * (sat**this%bc_eps(icell))
          this%wave_depth(2, icell) = this%ext_depth_uz(icell)
          theta2 = this%wave_theta(2, icell)
          flux2 = this%wave_flux(2, icell)
          flux1 = this%wave_flux(1, icell)
          theta1 = this%wave_theta(1, icell)
          this%wave_speed(2, icell) = &
            leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                      this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
          this%nwaves(icell) = this%nwaves(icell) + 1
          if (this%nwaves(icell) > this%nwaves_max) then
            !
            ! -- too many waves error
            ierr = 1
            goto 500
          end if
        end if
      end if
    else
      !
      ! -- extinction depth splits waves
      if (this%wave_depth(1, icell) - this%ext_depth_uz(icell) > DEM7) then
        jext = 2
        found = 0
        !
        ! -- locate extinction depth between waves
        do while (found == 0)
          diff = this%wave_depth(jext, icell) - this%ext_depth_uz(icell)
          if (diff > dzero) then
            jext = jext + 1
          else
            found = 1
          end if
        end do
        j = jext
        if (this%wave_theta(jext, icell) > theta_min) then
          !
          ! -- create a wave at extinction depth
          if (abs(diff) > DEM5) then
            if (this%nwaves(icell) + 1 > this%nwaves_max) then
              !
              ! -- too many waves error
              ierr = 1
              goto 500
            end if
            call this%shift_waves(icell, -1, &
                                  this%nwaves(icell) + 1, jext, -1)
            this%wave_depth(jext, icell) = this%ext_depth_uz(icell)
            this%nwaves(icell) = this%nwaves(icell) + 1
            if (this%nwaves(icell) > this%nwaves_max) then
              !
              ! -- too many waves error
              ierr = 1
              goto 500
            end if
          end if
          j = jext
        else
          jwet = this%nwaves(icell)
          j = jext + 1
          do while (j < this%nwaves(icell))
            if (this%wave_theta(j, icell) > theta_min) then
              jwet = j
              j = this%nwaves(icell) + 1
            end if
            j = j + 1
          end do
          jext = jwet
          j = jwet
        end if
      else
        j = 1
      end if
      !
      ! -- all waves above extinction depth
      do while (j <= this%nwaves(icell))
        if (ietflag == 2) then
          tho = this%wave_theta(j, icell)
          fktho = this%wave_flux(j, icell)
          hcap = this%caph(icell, tho)
          thetaout = this%rate_et_z(icell, factor, fktho, hcap)
        end if
        if (this%wave_theta(j, icell) > theta_min) then
          if (this%wave_theta(j, icell) - thetaout > &
              theta_min) then
            this%wave_theta(j, icell) = this%wave_theta(j, icell) - thetaout
          else if (this%wave_theta(j, icell) > theta_min) then
            this%wave_theta(j, icell) = theta_min
          end if
          if (j == 1) then
            this%wave_flux(j, icell) = &
              this%vks(icell) * &
              (((this%wave_theta(j, icell) - &
                 this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
          end if
          if (j > 1) then
            sat = (this%wave_theta(j - 1, icell) - this%theta_res(icell)) * &
                  dtheta_inv
            flux1 = this%vks(icell) * sat**this%bc_eps(icell)
            sat = (this%wave_theta(j, icell) - this%theta_res(icell)) * dtheta_inv
            flux2 = this%vks(icell) * sat**this%bc_eps(icell)
            this%wave_flux(j, icell) = flux2
            theta2 = this%wave_theta(j, icell)
            theta1 = this%wave_theta(j - 1, icell)
            this%wave_speed(j, icell) = &
              leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                        this%theta_res(icell), this%bc_eps(icell), &
                        this%vks(icell))
          end if
        end if
        j = j + 1
      end do
    end if
    !
    ! -- calculate aet
    j = 1
    do while (j <= this%nwaves(icell) - 1)
      if (abs(this%wave_theta(j, icell) - this%wave_theta(j + 1, &
                                                          icell)) < DEM6) then
        call this%shift_waves(icell, 1, j + 1, &
                              this%nwaves(icell) - 1, 1)
        j = j - 1
        this%nwaves(icell) = this%nwaves(icell) - 1
      end if
      j = j + 1
    end do
    depth = this%wave_depth(1, icell)
    fm = this%unsat_stor(icell, depth)
    this%et_uz(icell) = st - fm
    fm = this%et_uz(icell) / delt
    if (this%et_uz(icell) < dzero) then
      call this%load_waves(icell, this%etsav)
      this%nwaves(icell) = nwv
      this%et_uz(icell) = DZERO
    elseif (petsub - fm < -DEM15 .AND. ietflag == 2) then
      !
      ! -- aet greater than pet, reset and try again
      call this%load_waves(icell, this%etsav)
      this%nwaves(icell) = nwv
      this%et_uz(icell) = DZERO
    else
      itest = 1
    end if
    !
    ! -- end aet-pet loop for head dependent et
    fmp = fm
    if (k > 100) then
      itest = 1
    elseif (ietflag < 2) then
      fmp = petsub
      itest = 1
    end if
  end do
500 continue
  end procedure uzet

  !> @brief Calculate capillary pressure head from B-C equation
  !<
  module procedure caph
  ! -- local
  real(DP) :: lambda
  real(DP) :: star
  !
  caph = -DEM6
  star = (tho - this%theta_res(icell)) / &
         (this%theta_sat(icell) - this%theta_res(icell))
  if (star < DEM15) star = DEM15
  lambda = DTWO / (this%bc_eps(icell) - DTHREE)
  if (star > DEM15) then
    if (tho - this%theta_sat(icell) < DEM15) then
      caph = this%air_entry(icell) * star**(-DONE / lambda)
    else
      caph = DZERO
    end if
  end if
  end procedure caph

  !> @brief Calculate capillary pressure-based uz et
  !<
  module procedure rate_et_z
  rate_et_z = factor * fktho * (h - this%root_pot(icell))
  if (rate_et_z < DZERO) rate_et_z = DZERO
  end procedure rate_et_z

end submodule UzfCellGroupEt
