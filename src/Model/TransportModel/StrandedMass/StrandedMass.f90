!> @brief Stranded mass left behind by a falling water table
!!
!! Solute held by residual water and by the solid phase does not drain with the
!! water table. This module holds that mass out of the mobile system and returns
!! it when the cell rewets. Transfers are paired, so mass is conserved by
!! construction rather than by correction.
!<
module StrandedMassModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, DEM6, LENMEMPATH
  use MemoryManagerModule, only: mem_allocate, mem_deallocate
  use MemoryHelperModule, only: create_mem_path

  implicit none
  private
  public :: StrandedMassType
  public :: strand_rate, strand_rate_aqueous, strand_rate_sorbed
  public :: return_fraction, decay_amount

  !> @brief Stranded mass reservoirs for one transport domain
  !!
  !! The reservoir is split so that each part decays at the rate belonging to
  !! its phase. The owning package supplies its own properties; this type never
  !! reaches into a package.
  !<
  type :: StrandedMassType
    character(len=LENMEMPATH) :: memoryPath = '' !< memory path of the reservoir
    integer(I4B), pointer :: nodes => null() !< number of cells
    real(DP), dimension(:), pointer, contiguous :: stranded_aqueous => null() !< mass stranded from residual water
    real(DP), dimension(:), pointer, contiguous :: stranded_sorbed => null() !< mass stranded from the solid phase
    real(DP), dimension(:), pointer, contiguous :: ratestrand => null() !< mobile-side transfer rate, positive on return
    real(DP), dimension(:), pointer, contiguous :: ratedcystrand => null() !< decay rate of the aqueous reservoir
    real(DP), dimension(:), pointer, contiguous :: ratedcystrands => null() !< decay rate of the sorbed reservoir
  contains
    procedure :: init => strand_init
    procedure :: da => strand_da
    procedure :: total => strand_total
  end type StrandedMassType

contains

  !> @brief Allocate the reservoirs for a domain
  !<
  subroutine strand_init(this, name_model, subcomponent, nodes)
    class(StrandedMassType) :: this
    character(len=*), intent(in) :: name_model !< name of the transport model
    character(len=*), intent(in) :: subcomponent !< owning package, so paths do not collide
    integer(I4B), intent(in) :: nodes !< number of cells
    integer(I4B) :: n

    this%memoryPath = create_mem_path(name_model, subcomponent)

    call mem_allocate(this%nodes, 'NODES', this%memoryPath)
    this%nodes = nodes

    call mem_allocate(this%stranded_aqueous, nodes, 'STRANDED_AQUEOUS', &
                      this%memoryPath)
    call mem_allocate(this%stranded_sorbed, nodes, 'STRANDED_SORBED', &
                      this%memoryPath)
    call mem_allocate(this%ratestrand, nodes, 'RATESTRAND', this%memoryPath)
    call mem_allocate(this%ratedcystrand, nodes, 'RATEDCYSTRAND', &
                      this%memoryPath)
    call mem_allocate(this%ratedcystrands, nodes, 'RATEDCYSTRANDS', &
                      this%memoryPath)

    do n = 1, nodes
      this%stranded_aqueous(n) = DZERO
      this%stranded_sorbed(n) = DZERO
      this%ratestrand(n) = DZERO
      this%ratedcystrand(n) = DZERO
      this%ratedcystrands(n) = DZERO
    end do
  end subroutine strand_init

  !> @brief Deallocate the reservoirs
  !<
  subroutine strand_da(this)
    class(StrandedMassType) :: this

    call mem_deallocate(this%stranded_aqueous)
    call mem_deallocate(this%stranded_sorbed)
    call mem_deallocate(this%ratestrand)
    call mem_deallocate(this%ratedcystrand)
    call mem_deallocate(this%ratedcystrands)
    call mem_deallocate(this%nodes)
  end subroutine strand_da

  !> @brief Stranded mass held in cell n
  !<
  function strand_total(this, n) result(mass)
    class(StrandedMassType) :: this
    integer(I4B), intent(in) :: n !< cell number
    real(DP) :: mass

    mass = this%stranded_aqueous(n) + this%stranded_sorbed(n)
  end function strand_total

  !> @brief Rate at which solute in residual water is stranded by drainage
  !<
  pure function strand_rate_aqueous(ds, vcell, theta_r, conc, delt) result(rate)
    real(DP), intent(in) :: ds !< decrease in saturation over the step
    real(DP), intent(in) :: vcell !< cell volume
    real(DP), intent(in) :: theta_r !< residual water content
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: delt !< length of the time step
    real(DP) :: rate

    rate = ds * vcell * theta_r * conc / delt
  end function strand_rate_aqueous

  !> @brief Rate at which sorbed solute is stranded by drainage
  !<
  pure function strand_rate_sorbed(ds, vcell, volfracm, rhob, sval, delt) &
    result(rate)
    real(DP), intent(in) :: ds !< decrease in saturation over the step
    real(DP), intent(in) :: vcell !< cell volume
    real(DP), intent(in) :: volfracm !< volume fraction of the domain
    real(DP), intent(in) :: rhob !< bulk density
    real(DP), intent(in) :: sval !< sorbed concentration from the isotherm
    real(DP), intent(in) :: delt !< length of the time step
    real(DP) :: rate

    rate = ds * vcell * volfracm * rhob * sval / delt
  end function strand_rate_sorbed

  !> @brief Total rate at which mass is stranded by drainage
  !<
  pure function strand_rate(ds, vcell, theta_r, conc, volfracm, rhob, sval, &
                            delt) result(rate)
    real(DP), intent(in) :: ds !< decrease in saturation over the step
    real(DP), intent(in) :: vcell !< cell volume
    real(DP), intent(in) :: theta_r !< residual water content
    real(DP), intent(in) :: conc !< solute concentration
    real(DP), intent(in) :: volfracm !< volume fraction of the domain
    real(DP), intent(in) :: rhob !< bulk density
    real(DP), intent(in) :: sval !< sorbed concentration from the isotherm
    real(DP), intent(in) :: delt !< length of the time step
    real(DP) :: rate

    rate = strand_rate_aqueous(ds, vcell, theta_r, conc, delt) + &
           strand_rate_sorbed(ds, vcell, volfracm, rhob, sval, delt)
  end function strand_rate

  !> @brief Share of the reservoir returned by rewetting
  !!
  !! The stranded mass is taken to be spread uniformly through the unsaturated
  !! part of the cell, so the newly rewetted share returns.
  !<
  pure function return_fraction(dw, sat_old) result(f)
    real(DP), intent(in) :: dw !< increase in saturation over the step
    real(DP), intent(in) :: sat_old !< saturation at the end of the last step
    real(DP) :: f
    real(DP) :: unsat

    unsat = DONE - sat_old
    if (unsat < DEM6 .or. dw <= DZERO) then
      f = DZERO
    else
      f = dw / unsat
      if (f > DONE) f = DONE
    end if
  end function return_fraction

  !> @brief Mass lost from a reservoir to first-order decay over a step
  !!
  !! A negative rate is production, following the sign convention of the decay
  !! rates of the mobile domain, and returns a negative amount. Production is
  !! proportional to the mass held, so an empty reservoir stays empty.
  !<
  pure function decay_amount(mass, lambda, delt) result(amount)
    real(DP), intent(in) :: mass !< mass held in the reservoir
    real(DP), intent(in) :: lambda !< first-order decay rate
    real(DP), intent(in) :: delt !< length of the time step
    real(DP) :: amount

    if (lambda == DZERO .or. mass <= DZERO) then
      amount = DZERO
    else
      amount = mass * (DONE - exp(-lambda * delt))
    end if
  end function decay_amount

end module StrandedMassModule
