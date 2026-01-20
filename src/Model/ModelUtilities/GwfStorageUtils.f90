!> @brief This module contains stateless storage subroutines and functions
!!
!! This module contains the functions to calculate the specific
!! storage (SC1) and specific yield (SC2) capacities that are used in
!! the storage (STO) package. It also contains subroutines to calculate
!! the amat and rhs terms for specific storage and specific yield.
!! This module does not depend on the STO package.
!!
!<
module GwfStorageUtilsModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DHALF, DONE

  implicit none
  private
  public :: SsTerms
  public :: SyTerms
  public :: SsCapacity
  public :: SyCapacity

contains

  !> @brief Calculate the specific storage terms
  !!
  !! Subroutine to calculate the specific storage terms for a cell using
  !! the cell geometry, current and previous specific storage capacity,
  !! current and previous cell saturation, and the current and previous head.
  !! Subroutine can optionally return the flow rate from specific storage.
  !!
  !<
  pure subroutine SsTerms(iconvert, iorig_ss, iconf_ss, top, bot, &
                          rho1, rho1old, snnew, snold, hnew, hold, &
                          aterm, rhsterm, rate)
    ! -- dummy variables
    integer(I4B), intent(in) :: iconvert !< flag indicating if cell is convertible
    integer(I4B), intent(in) :: iorig_ss !< flag indicating if the original MODFLOW 6 specific storage formulation is being used
    integer(I4B), intent(in) :: iconf_ss !< flag indicating if specific storage only applies under confined conditions
    real(DP), intent(in) :: top !< top of cell
    real(DP), intent(in) :: bot !< bottom of cell
    real(DP), intent(in) :: rho1 !< current specific storage capacity
    real(DP), intent(in) :: rho1old !< previous specific storage capacity
    real(DP), intent(in) :: snnew !< current cell saturation
    real(DP), intent(in) :: snold !< previous cell saturation
    real(DP), intent(in) :: hnew !< current head
    real(DP), intent(in) :: hold !< previous head
    real(DP), intent(inout) :: aterm !< coefficient matrix term
    real(DP), intent(inout) :: rhsterm !< right-hand side term
    real(DP), intent(inout), optional :: rate !< calculated specific storage rate
    ! -- local variables
    real(DP) :: tthk
    real(DP) :: zold
    real(DP) :: znew
    !
    ! -- initialize terms
    aterm = -rho1 * snnew
    rhsterm = DZERO
    !
    ! -- calculate specific storage terms
    if (iconvert /= 0) then
      if (iorig_ss == 0) then
        if (iconf_ss == 0) then
          tthk = top - bot
          zold = bot + DHALF * tthk * snold
          znew = bot + DHALF * tthk * snnew
          rhsterm = -rho1old * snold * (hold - zold) - rho1 * snnew * znew
        else
          if (snold == DONE) then
            rhsterm = rhsterm - rho1old * (hold - top)
          end if
          if (snnew == DONE) then
            rhsterm = rhsterm - rho1 * top
          else
            aterm = DZERO
          end if
        end if
      else
        rhsterm = -rho1old * snold * hold
      end if
    else
      rhsterm = -rho1old * snold * hold
    end if
    !
    ! -- calculate rate
    if (present(rate)) then
      rate = aterm * hnew - rhsterm
    end if
  end subroutine SsTerms

  !> @brief Calculate the specific yield storage terms
  !!
  !! Subroutine to calculate the specific yield storage terms for a cell
  !! using the cell geometry, current and previous specific yield storage
  !! capacity, and the current and previous cell saturation. Subroutine
  !! can optionally return the flow rate from specific yield.
  !!
  !<
  pure subroutine SyTerms(top, bot, rho2, rho2old, snnew, snold, &
                          aterm, rhsterm, rate)
    ! -- dummy variables
    real(DP), intent(in) :: top !< top of cell
    real(DP), intent(in) :: bot !< bottom of cell
    real(DP), intent(in) :: rho2 !< current specific yield storage capacity
    real(DP), intent(in) :: rho2old !< previous specific yield storage capacity
    real(DP), intent(in) :: snnew !< current cell saturation
    real(DP), intent(in) :: snold !< previous cell saturation
    real(DP), intent(inout) :: aterm !< coefficient matrix term
    real(DP), intent(inout) :: rhsterm !< right-hand side term
    real(DP), intent(inout), optional :: rate !< calculated specific yield rate
    ! -- local variables
    real(DP) :: tthk
    !
    ! -- initialize terms
    aterm = DZERO
    tthk = top - bot
    !
    ! -- calculate specific yield storage terms
    if (snnew < DONE) then
      if (snnew > DZERO) then
        aterm = -rho2
        rhsterm = -rho2old * tthk * snold - rho2 * bot
      else
        rhsterm = tthk * (DZERO - rho2old * snold)
      end if
      ! -- known flow from specific yield
    else
      rhsterm = tthk * (rho2 * snnew - rho2old * snold)
    end if
    !
    ! -- calculate rate
    if (present(rate)) then
      rate = rho2old * tthk * snold - rho2 * tthk * snnew
    end if
  end subroutine SyTerms

  !> @brief Calculate the specific storage capacity
  !!
  !! Function to calculate the specific storage capacity using
  !! the cell geometry and the specific storage or storage coefficient.
  !!
  !! @return      sc1               specific storage capacity
  !<
  pure function SsCapacity(istor_coef, top, bot, area, ss) result(sc1)
    ! -- dummy variables
    integer(I4B), intent(in) :: istor_coef !< flag indicating if ss is the storage coefficient
    real(DP), intent(in) :: top !< top of cell
    real(DP), intent(in) :: bot !< bottom of cell
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: ss !< specific storage or storage coefficient
    ! -- local variables
    real(DP) :: sc1
    real(DP) :: thick
    ! -- calculate specific storage capacity
    if (istor_coef == 0) then
      thick = top - bot
    else
      thick = DONE
    end if
    sc1 = ss * thick * area
  end function SsCapacity

  !> @brief Calculate the specific yield capacity
  !!
  !! Function to calculate the specific yield capacity using
  !! the cell area and the specific yield.
  !!
  !! @return      sc2               specific yield capacity
  !<
  pure function SyCapacity(area, sy) result(sc2)
    ! -- dummy variables
    real(DP), intent(in) :: area !< horizontal cell area
    real(DP), intent(in) :: sy !< specific yield
    ! -- local variables
    real(DP) :: sc2
    ! -- calculate specific yield capacity
    sc2 = sy * area
  end function SyCapacity

end module GwfStorageUtilsModule
