  !> @brief Kinematic routing equation residual
    !!
    !! Method to calculate the kinematic-wave routing
    !! residual.
    !!
  !<
  function kinematic_residual(qa, qb, qc, qd, &
                              aa, ab, ac, ad, &
                              qsrc, length, weight, delt, &
                              courant)
    use KindModule, only: DP
    use ConstantsModule, only: DZERO, DHALF, DONE
    ! -- return variable
    real(DP) :: kinematic_residual !< kinematic-wave residual
    ! -- dummy variables
    real(DP), intent(in) :: qa !< upstream flow at previous time
    real(DP), intent(in) :: qb !< downstream flow at previous time
    real(DP), intent(in) :: qc !< upstream flow
    real(DP), intent(in) :: qd !< downstream flow
    real(DP), intent(in) :: aa !< upstream area at previous time
    real(DP), intent(in) :: ab !< downstream area at previous time
    real(DP), intent(in) :: ac !< upstream area
    real(DP), intent(in) :: ad !< downstream area
    real(DP), intent(in) :: qsrc !< lateral flow term (L3T-1L-1)
    real(DP), intent(in) :: length !< reach length
    real(DP), intent(in) :: weight !< temporal weight
    real(DP), intent(in) :: delt !< time step length
    real(DP), intent(in) :: courant !< courant number
    ! --local variables
    real(DP) :: f11
    real(DP) :: f12

    if (courant > DONE) then
      f11 = (qd - qc) / length
      f12 = (ac - aa) / delt
    else
      f11 = (weight * (qd - qc) + (DONE - weight) * (qb - qa)) / length
      f12 = DHALF * ((ad - ab) + (ac - aa)) / delt
    end if
    kinematic_residual = f11 + f12 - qsrc

  end function kinematic_residual

  !> @brief Kinematic routing equation storage term
    !!
    !! Method to calculate the kinematic-wave routing
    !! storage term.
    !!
  !<
  function kinematic_storage(aa, ab, ac, ad, length, delt, courant)
    use KindModule, only: DP
    use ConstantsModule, only: DHALF

    real(DP) :: kinematic_storage !< kinematic-wave storage change

    real(DP), intent(in) :: aa !< upstream area at previous time
    real(DP), intent(in) :: ab !< downstream area at previous time
    real(DP), intent(in) :: ac !< upstream area
    real(DP), intent(in) :: ad !< downstream area
    real(DP), intent(in) :: length !< reach length
    real(DP), intent(in) :: delt !< time step length
    real(DP), intent(in) :: courant !< courant number

    if (courant > DONE) then
      kinematic_storage = (aa - ac) * length / delt
    else
      kinematic_storage = DHALF * ((ac + ad) - (aa + ab)) * length / delt
    end if

  end function kinematic_storage
