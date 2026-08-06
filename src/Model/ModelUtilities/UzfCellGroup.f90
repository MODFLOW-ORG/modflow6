!> @brief Unsaturated zone flow by the method of characteristics
!!
!! Water moves through the unsaturated zone as a train of kinematic waves, after
!! Smith (1983) and the MODFLOW-2005 UZF1 Package (Niswonger and others, 2006).
!! A wave is a step change in water content moving down at speed dflux/dtheta,
!! with flux given by the Brooks-Corey relation
!! flux = vks * ((theta - theta_res) / (theta_sat - theta_res))**bc_eps.
!!
!! The waves of one cell are held in columns of wave_depth, wave_theta,
!! wave_flux and wave_speed, packed into positions 1 to nwaves(icell). The
!! train is ordered so that wave 1 is deepest, depth decreases with index, and
!! the wave at nwaves is the one at the land surface. Increasing infiltration
!! adds a lead wave at the surface; decreasing infiltration adds a set of
!! ntrailwaves trailing waves. Waves are routed until one reaches the bottom of
!! the unsaturated zone, where it becomes recharge, or overtakes the wave ahead
!! of it, where the two merge. Both events remove a wave and compact the train.
!<
module UzfCellGroupModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DEM30, DEM20, DEM15, DEM14, DEM12, DEM10, &
                             DEM9, DEM7, DEM6, DEM5, DEM4, DEM3, DHALF, DONE, &
                             DTWO, DTHREE, DEP20
  use SmoothingModule
  use TdisModule, only: ITMUNI, delt, kper
  use UzfETUtilModule, only: etfunc_lin, etfunc_nlin

  implicit none
  private
  public :: UzfCellGroupType

  !> @brief Snapshot of one cell's wave train
  !!
  !! Holds the waves of a single cell so they can be restored after a trial
  !! solution. Two are kept: one for the outer-iteration reset in solve() and
  !! one for the aet-to-pet retry loop in uzet().
  !<
  type :: UzfWaveStoreType
    real(DP), pointer, dimension(:), contiguous :: depth => null()
    real(DP), pointer, dimension(:), contiguous :: theta => null()
    real(DP), pointer, dimension(:), contiguous :: flux => null()
    real(DP), pointer, dimension(:), contiguous :: speed => null()
    integer(I4B), pointer :: nwaves => null()
  end type UzfWaveStoreType

  type :: UzfCellGroupType
    !
    ! -- wave train; dimensioned (nwaves_max, ncells), with the waves of cell
    !    icell occupying positions 1 to nwaves(icell)
    real(DP), dimension(:, :), pointer, contiguous :: wave_depth => null() !< depth of the wave front below celtop
    real(DP), dimension(:, :), pointer, contiguous :: wave_theta => null() !< water content behind the wave front
    real(DP), dimension(:, :), pointer, contiguous :: wave_flux => null() !< unsaturated flux behind the wave front
    real(DP), dimension(:, :), pointer, contiguous :: wave_speed => null() !< wave speed, dflux/dtheta
    integer(I4B), pointer, dimension(:), contiguous :: nwaves => null() !< number of active waves in each cell
    integer(I4B), pointer :: nwaves_max => null() !< wave capacity, NTRAILWAVES * NWAVESETS
    integer(I4B), pointer :: ntrailwaves => null() !< trailing waves created when infiltration drops
    !
    ! -- unsaturated zone properties
    real(DP), pointer, dimension(:), contiguous :: theta_res => null() !< residual water content
    real(DP), pointer, dimension(:), contiguous :: theta_sat => null() !< saturated water content
    real(DP), pointer, dimension(:), contiguous :: theta_init => null() !< initial water content
    real(DP), pointer, dimension(:), contiguous :: bc_eps => null() !< Brooks-Corey exponent
    real(DP), pointer, dimension(:), contiguous :: vks => null() !< vertical saturated hydraulic conductivity
    !
    ! -- geometry
    real(DP), pointer, dimension(:), contiguous :: celtop => null() !< top of the unsaturated zone in this cell
    real(DP), pointer, dimension(:), contiguous :: celbot => null() !< bottom of this cell
    real(DP), pointer, dimension(:), contiguous :: landtop => null() !< land surface above this cell
    real(DP), pointer, dimension(:), contiguous :: surfdep => null() !< undulation depth of the land surface
    real(DP), pointer, dimension(:), contiguous :: cellarea => null() !< area of the host gwf cell
    real(DP), pointer, dimension(:), contiguous :: uzfarea => null() !< area of this uzf cell
    integer(I4B), pointer, dimension(:), contiguous :: landflag => null() !< cell is at land surface (1) or not (0)
    integer(I4B), pointer, dimension(:), contiguous :: cell_below => null() !< uzf cell directly below, 0 if none
    !
    ! -- water table
    real(DP), pointer, dimension(:), contiguous :: water_table => null() !< water table elevation this iteration
    real(DP), pointer, dimension(:), contiguous :: water_table_old => null() !< water table elevation last time step
    !
    ! -- infiltration and discharge
    real(DP), pointer, dimension(:), contiguous :: finf_spec => null() !< specified infiltration rate
    real(DP), pointer, dimension(:), contiguous :: finf => null() !< infiltration rate applied to this cell
    real(DP), pointer, dimension(:), contiguous :: finf_rej => null() !< infiltration rejected at land surface
    real(DP), pointer, dimension(:), contiguous :: surf_infil => null() !< infiltration entering the unsaturated zone
    real(DP), pointer, dimension(:), contiguous :: surf_infil_below => null() !< infiltration passed to cell_below
    real(DP), pointer, dimension(:), contiguous :: surf_seep => null() !< groundwater discharge to land surface
    real(DP), pointer, dimension(:), contiguous :: flux_to_wt => null() !< water reaching the water table this time step
    !
    ! -- evapotranspiration
    real(DP), pointer, dimension(:), contiguous :: pet => null() !< potential et available to this cell
    real(DP), pointer, dimension(:), contiguous :: pet_max => null() !< potential et at land surface
    real(DP), pointer, dimension(:), contiguous :: gw_pet => null() !< potential et left for the saturated zone
    real(DP), pointer, dimension(:), contiguous :: gwet => null() !< et taken from the saturated zone
    real(DP), pointer, dimension(:), contiguous :: et_uz => null() !< et taken from the unsaturated zone
    real(DP), pointer, dimension(:), contiguous :: ext_depth => null() !< extinction depth below land surface
    real(DP), pointer, dimension(:), contiguous :: ext_depth_uz => null() !< part of the extinction depth in this cell
    real(DP), pointer, dimension(:), contiguous :: theta_ext => null() !< extinction water content
    real(DP), pointer, dimension(:), contiguous :: air_entry => null() !< air entry potential
    real(DP), pointer, dimension(:), contiguous :: root_pot => null() !< root potential
    real(DP), pointer, dimension(:), contiguous :: root_act => null() !< root activity function
    !
    ! -- wave snapshots, one cell each
    type(UzfWaveStoreType) :: wavsav !< waves saved across an outer iteration
    type(UzfWaveStoreType) :: etsav !< waves saved across the uzet retry loop
    !
    ! -- work arrays for leadwav, sized for the wave capacity of one cell
    real(DP), pointer, dimension(:), contiguous :: overtake_time => null() !< time until each wave overtakes the one ahead
    integer(I4B), pointer, dimension(:), contiguous :: wave_merges => null() !< wave merges with the one ahead this substep

  contains

    procedure :: init
    procedure :: setdata
    procedure :: sethead
    procedure :: setdatauzfarea
    procedure :: setdatafinf
    procedure :: setdataet
    procedure :: setdataetwc
    procedure :: setdataetha
    procedure :: setwaves
    procedure :: shift_waves
    procedure :: store_waves
    procedure :: load_waves
    procedure :: routewaves
    procedure :: uzflow
    procedure :: addrech
    procedure :: trailwav
    procedure :: leadwav
    procedure :: advance
    procedure :: solve
    procedure :: unsat_stor
    procedure :: update_wav
    procedure :: simgwet
    procedure :: caph
    procedure :: rate_et_z
    procedure :: uzet
    procedure :: uz_rise
    procedure :: rejfinf
    procedure :: gwseep
    procedure :: setbelowpet
    procedure :: setgwpet
    procedure :: dealloc
    procedure :: get_water_content_at_depth
    procedure :: get_wcnew
  end type UzfCellGroupType

  interface

    module subroutine setwaves(this, icell)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
    end subroutine setwaves

    module subroutine routewaves(this, totfluxtot, delt, ietflag, icell, ierr)
      class(UzfCellGroupType) :: this
      real(DP), intent(inout) :: totfluxtot
      real(DP), intent(in) :: delt
      integer(I4B), intent(in) :: ietflag
      integer(I4B), intent(in) :: icell
      integer(I4B), intent(inout) :: ierr
    end subroutine routewaves

    module subroutine shift_waves(this, icell, shft, strt, stp, cntr)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      integer(I4B), intent(in) :: shft
      integer(I4B), intent(in) :: strt
      integer(I4B), intent(in) :: stp
      integer(I4B), intent(in) :: cntr
    end subroutine shift_waves

    module subroutine store_waves(this, icell, store)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      type(UzfWaveStoreType), intent(inout) :: store
    end subroutine store_waves

    module subroutine load_waves(this, icell, store)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      type(UzfWaveStoreType), intent(in) :: store
    end subroutine load_waves

    module subroutine uzflow(this, thick, thickold, delt, ietflag, icell, ierr)
      class(UzfCellGroupType) :: this
      real(DP), intent(inout) :: thickold
      real(DP), intent(inout) :: thick
      real(DP), intent(in) :: delt
      integer(I4B), intent(in) :: ietflag
      integer(I4B), intent(in) :: icell
      integer(I4B), intent(inout) :: ierr
    end subroutine uzflow

    module subroutine factors(eps_thick, eps_flux)
      real(DP), intent(out) :: eps_thick
      real(DP), intent(out) :: eps_flux
    end subroutine factors

    module subroutine trailwav(this, icell, ierr)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      integer(I4B), intent(inout) :: ierr
    end subroutine trailwav

  module subroutine leadwav(this, time, single_wave, trail_added, thetab, fluxb, &
                              dflux_surf, eps_flux, delt, icell)
      class(UzfCellGroupType) :: this
      real(DP), intent(inout) :: thetab
      real(DP), intent(inout) :: fluxb
      real(DP), intent(in) :: eps_flux
      real(DP), intent(inout) :: time
      integer(I4B), intent(inout) :: single_wave
      integer(I4B), intent(inout) :: trail_added
      real(DP), intent(inout) :: dflux_surf
      real(DP), intent(in) :: delt
      integer(I4B), intent(in) :: icell
    end subroutine leadwav

    module function unsat_stor(this, icell, d1)
      real(DP) :: unsat_stor
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      real(DP), intent(inout) :: d1
    end function unsat_stor

    module subroutine update_wav(this, icell, delt, iss, itest)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      integer(I4B), intent(in) :: itest
      integer(I4B), intent(in) :: iss
      real(DP), intent(in) :: delt
    end subroutine update_wav

    module subroutine uz_rise(this, icell, totfluxtot)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      real(DP), intent(inout) :: totfluxtot
    end subroutine uz_rise

    module function get_water_content_at_depth(this, icell, depth) result(theta_at_depth)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell !< uzf cell containing depth
      real(DP), intent(in) :: depth !< depth within the cell
      real(DP) :: theta_at_depth
    end function get_water_content_at_depth

    module function get_wcnew(this, icell) result(watercontent)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell !< uzf cell containing depth
      real(DP) :: watercontent
    end function get_wcnew

    module function leadspeed(theta1, theta2, flux1, flux2, thts, thtr, eps, vks)
      real(DP) :: leadspeed
      real(DP), intent(in) :: theta1
      real(DP), intent(in) :: theta2
      real(DP), intent(in) :: flux1
      real(DP), intent(inout) :: flux2
      real(DP), intent(in) :: thts
      real(DP), intent(in) :: thtr
      real(DP), intent(in) :: eps
      real(DP), intent(in) :: vks
    end function leadspeed

    module subroutine setgwpet(this, icell)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
    end subroutine setgwpet

    module subroutine setbelowpet(this, icell, jbelow)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      integer(I4B), intent(in) :: jbelow
    end subroutine setbelowpet

    module subroutine simgwet(this, igwetflag, icell, hgwf, trhs, thcof, det)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: igwetflag
      integer(I4B), intent(in) :: icell
      real(DP), intent(in) :: hgwf
      real(DP), intent(inout) :: trhs
      real(DP), intent(inout) :: thcof
      real(DP), intent(inout) :: det
    end subroutine simgwet

    module subroutine uzet(this, icell, delt, ietflag, ierr)
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      real(DP), intent(in) :: delt
      integer(I4B), intent(in) :: ietflag
      integer(I4B), intent(inout) :: ierr
    end subroutine uzet

    module function caph(this, icell, tho)
      real(DP) :: caph
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      real(DP), intent(in) :: tho
    end function caph

    module function rate_et_z(this, icell, factor, fktho, h)
      real(DP) :: rate_et_z
      class(UzfCellGroupType) :: this
      integer(I4B), intent(in) :: icell
      real(DP), intent(in) :: factor, fktho, h
    end function rate_et_z

  end interface

contains

  !> @brief Allocate and set uzf object variables
  !<
  subroutine init(this, ncells, nwav, memory_path)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: ncells
    integer(I4B), intent(in) :: nwav
    character(len=*), intent(in) :: memory_path
    ! -- local
    integer(I4B) :: icell
    integer(I4B) :: j
    !
    ! -- wave state, one column per cell
    call mem_allocate(this%wave_depth, nwav, ncells, 'WAVE_DEPTH', memory_path)
    call mem_allocate(this%wave_theta, nwav, ncells, 'WAVE_THETA', memory_path)
    call mem_allocate(this%wave_flux, nwav, ncells, 'WAVE_FLUX', memory_path)
    call mem_allocate(this%wave_speed, nwav, ncells, 'WAVE_SPEED', memory_path)
    call mem_allocate(this%nwaves, ncells, 'NWAVES', memory_path)
    call mem_allocate(this%nwaves_max, 'NWAVES_MAX', memory_path)
    call mem_allocate(this%ntrailwaves, 'NTRAILWAVES', memory_path)
    !
    ! -- cell properties and state
    call mem_allocate(this%theta_res, ncells, 'THETA_RES', memory_path)
    call mem_allocate(this%theta_sat, ncells, 'THETA_SAT', memory_path)
    call mem_allocate(this%theta_init, ncells, 'THETA_INIT', memory_path)
    call mem_allocate(this%bc_eps, ncells, 'BC_EPS', memory_path)
    call mem_allocate(this%air_entry, ncells, 'AIR_ENTRY', memory_path)
    call mem_allocate(this%root_pot, ncells, 'ROOT_POT', memory_path)
    call mem_allocate(this%root_act, ncells, 'ROOT_ACT', memory_path)
    call mem_allocate(this%theta_ext, ncells, 'THETA_EXT', memory_path)
    call mem_allocate(this%et_uz, ncells, 'ET_UZ', memory_path)
    call mem_allocate(this%flux_to_wt, ncells, 'FLUX_TO_WT', memory_path)
    call mem_allocate(this%finf_spec, ncells, 'FINF_SPEC', memory_path)
    call mem_allocate(this%finf, ncells, 'FINF', memory_path)
    call mem_allocate(this%finf_rej, ncells, 'FINF_REJ', memory_path)
    call mem_allocate(this%gwet, ncells, 'GWET', memory_path)
    call mem_allocate(this%uzfarea, ncells, 'UZFAREA', memory_path)
    call mem_allocate(this%cellarea, ncells, 'CELLAREA', memory_path)
    call mem_allocate(this%celtop, ncells, 'CELTOP', memory_path)
    call mem_allocate(this%celbot, ncells, 'CELBOT', memory_path)
    call mem_allocate(this%landtop, ncells, 'LANDTOP', memory_path)
    call mem_allocate(this%water_table, ncells, 'WATER_TABLE', memory_path)
   call mem_allocate(this%water_table_old, ncells, 'WATER_TABLE_OLD', memory_path)
    call mem_allocate(this%surfdep, ncells, 'SURFDEP', memory_path)
    call mem_allocate(this%vks, ncells, 'VKS', memory_path)
    call mem_allocate(this%surf_infil, ncells, 'SURF_INFIL', memory_path)
 call mem_allocate(this%surf_infil_below, ncells, 'SURF_INFIL_BELOW', memory_path)
    call mem_allocate(this%surf_seep, ncells, 'SURF_SEEP', memory_path)
    call mem_allocate(this%gw_pet, ncells, 'GW_PET', memory_path)
    call mem_allocate(this%pet, ncells, 'PET', memory_path)
    call mem_allocate(this%pet_max, ncells, 'PET_MAX', memory_path)
    call mem_allocate(this%ext_depth, ncells, 'EXT_DEPTH', memory_path)
    call mem_allocate(this%ext_depth_uz, ncells, 'EXT_DEPTH_UZ', memory_path)
    call mem_allocate(this%landflag, ncells, 'LANDFLAG', memory_path)
    call mem_allocate(this%cell_below, ncells, 'CELL_BELOW', memory_path)
    !
    ! -- wave snapshots and leadwav work arrays, one cell wide
    call mem_allocate(this%wavsav%depth, nwav, 'WSAV_DEPTH', memory_path)
    call mem_allocate(this%wavsav%theta, nwav, 'WSAV_THETA', memory_path)
    call mem_allocate(this%wavsav%flux, nwav, 'WSAV_FLUX', memory_path)
    call mem_allocate(this%wavsav%speed, nwav, 'WSAV_SPEED', memory_path)
    call mem_allocate(this%wavsav%nwaves, 'WSAV_NWAVES', memory_path)
    call mem_allocate(this%etsav%depth, nwav, 'ETSAV_DEPTH', memory_path)
    call mem_allocate(this%etsav%theta, nwav, 'ETSAV_THETA', memory_path)
    call mem_allocate(this%etsav%flux, nwav, 'ETSAV_FLUX', memory_path)
    call mem_allocate(this%etsav%speed, nwav, 'ETSAV_SPEED', memory_path)
    call mem_allocate(this%etsav%nwaves, 'ETSAV_NWAVES', memory_path)
    call mem_allocate(this%overtake_time, nwav, 'OVERTAKE_TIME', memory_path)
    call mem_allocate(this%wave_merges, nwav, 'WAVE_MERGES', memory_path)
    !
    ! -- wave capacity and trailing wave count are the same for every cell
    this%nwaves_max = nwav
    this%ntrailwaves = 0
    this%wavsav%nwaves = 0
    this%etsav%nwaves = 0
    do j = 1, nwav
      this%wavsav%depth(j) = DZERO
      this%wavsav%theta(j) = DZERO
      this%wavsav%flux(j) = DZERO
      this%wavsav%speed(j) = DZERO
      this%etsav%depth(j) = DZERO
      this%etsav%theta(j) = DZERO
      this%etsav%flux(j) = DZERO
      this%etsav%speed(j) = DZERO
      this%overtake_time(j) = DZERO
      this%wave_merges(j) = 0
    end do
    do icell = 1, ncells
      do j = 1, nwav
        this%wave_depth(j, icell) = DZERO
        this%wave_theta(j, icell) = DZERO
        this%wave_flux(j, icell) = DZERO
        this%wave_speed(j, icell) = DZERO
      end do
      this%nwaves(icell) = 1
      this%theta_res(icell) = DZERO
      this%theta_sat(icell) = DZERO
      this%theta_init(icell) = DZERO
      this%bc_eps(icell) = DZERO
      this%air_entry(icell) = DZERO
      this%root_pot(icell) = DZERO
      this%root_act(icell) = DZERO
      this%theta_ext(icell) = DZERO
      this%et_uz(icell) = DZERO
      this%flux_to_wt(icell) = DZERO
      this%finf_spec(icell) = DZERO
      this%finf(icell) = DZERO
      this%finf_rej(icell) = DZERO
      this%gwet(icell) = DZERO
      this%uzfarea(icell) = DZERO
      this%cellarea(icell) = DZERO
      this%celtop(icell) = DZERO
      this%celbot(icell) = DZERO
      this%landtop(icell) = DZERO
      this%water_table(icell) = DZERO
      this%water_table_old(icell) = DZERO
      this%surfdep(icell) = DZERO
      this%vks(icell) = DZERO
      this%surf_infil(icell) = DZERO
      this%surf_infil_below(icell) = DZERO
      this%surf_seep(icell) = DZERO
      this%gw_pet(icell) = DZERO
      this%pet(icell) = DZERO
      this%pet_max(icell) = DZERO
      this%ext_depth(icell) = DZERO
      this%ext_depth_uz(icell) = DZERO
      this%landflag(icell) = 0
      this%cell_below(icell) = 0
    end do
  end subroutine init

  !> @brief Deallocate uzf object variables
  !<
  subroutine dealloc(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(UzfCellGroupType) :: this
    !
    call mem_deallocate(this%wave_depth)
    call mem_deallocate(this%wave_theta)
    call mem_deallocate(this%wave_flux)
    call mem_deallocate(this%wave_speed)
    call mem_deallocate(this%nwaves)
    call mem_deallocate(this%nwaves_max)
    call mem_deallocate(this%ntrailwaves)
    call mem_deallocate(this%theta_res)
    call mem_deallocate(this%theta_sat)
    call mem_deallocate(this%theta_init)
    call mem_deallocate(this%bc_eps)
    call mem_deallocate(this%air_entry)
    call mem_deallocate(this%root_pot)
    call mem_deallocate(this%root_act)
    call mem_deallocate(this%theta_ext)
    call mem_deallocate(this%et_uz)
    call mem_deallocate(this%flux_to_wt)
    call mem_deallocate(this%finf_spec)
    call mem_deallocate(this%finf)
    call mem_deallocate(this%finf_rej)
    call mem_deallocate(this%gwet)
    call mem_deallocate(this%uzfarea)
    call mem_deallocate(this%cellarea)
    call mem_deallocate(this%celtop)
    call mem_deallocate(this%celbot)
    call mem_deallocate(this%landtop)
    call mem_deallocate(this%water_table)
    call mem_deallocate(this%water_table_old)
    call mem_deallocate(this%surfdep)
    call mem_deallocate(this%vks)
    call mem_deallocate(this%surf_infil)
    call mem_deallocate(this%surf_infil_below)
    call mem_deallocate(this%surf_seep)
    call mem_deallocate(this%gw_pet)
    call mem_deallocate(this%pet)
    call mem_deallocate(this%pet_max)
    call mem_deallocate(this%ext_depth)
    call mem_deallocate(this%ext_depth_uz)
    call mem_deallocate(this%landflag)
    call mem_deallocate(this%cell_below)
    call mem_deallocate(this%wavsav%depth)
    call mem_deallocate(this%wavsav%theta)
    call mem_deallocate(this%wavsav%flux)
    call mem_deallocate(this%wavsav%speed)
    call mem_deallocate(this%wavsav%nwaves)
    call mem_deallocate(this%etsav%depth)
    call mem_deallocate(this%etsav%theta)
    call mem_deallocate(this%etsav%flux)
    call mem_deallocate(this%etsav%speed)
    call mem_deallocate(this%etsav%nwaves)
    call mem_deallocate(this%overtake_time)
    call mem_deallocate(this%wave_merges)
  end subroutine dealloc

  !> @brief Set uzf object material properties
  !<
  subroutine setdata(this, icell, area, top, bot, surfdep, vks, thtr, thts, &
                     thti, eps, ntrail, landflag, ivertcon)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: area
    real(DP), intent(in) :: top
    real(DP), intent(in) :: bot
    real(DP), intent(in) :: surfdep
    real(DP), intent(in) :: vks
    real(DP), intent(in) :: thtr
    real(DP), intent(in) :: thts
    real(DP), intent(in) :: thti
    real(DP), intent(in) :: eps
    integer(I4B), intent(in) :: ntrail
    integer(I4B), intent(in) :: landflag
    integer(I4B), intent(in) :: ivertcon
    !
    ! -- set the values for uzf cell icell
    this%landflag(icell) = landflag
    this%cell_below(icell) = ivertcon
    this%surfdep(icell) = surfdep
    this%uzfarea(icell) = area
    this%cellarea(icell) = area
    if (this%landflag(icell) == 1) then
      this%celtop(icell) = top - DHALF * this%surfdep(icell)
    else
      this%celtop(icell) = top
    end if
    this%celbot(icell) = bot
    this%vks(icell) = vks
    this%theta_res(icell) = thtr
    this%theta_sat(icell) = thts
    this%theta_init(icell) = thti
    this%bc_eps(icell) = eps
    this%ntrailwaves = ntrail
    this%pet(icell) = DZERO
    this%ext_depth(icell) = DZERO
    this%theta_ext(icell) = DZERO
    this%air_entry(icell) = DZERO
    this%root_pot(icell) = DZERO
  end subroutine setdata

  !> @brief Set initial head for uzf object
  !<
  subroutine sethead(this, icell, hgwf)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: hgwf
    !
    ! -- set initial head
    this%water_table(icell) = this%celbot(icell)
    if (hgwf > this%celbot(icell)) this%water_table(icell) = hgwf
    if (this%water_table(icell) > this%celtop(icell)) &
      this%water_table(icell) = this%celtop(icell)
    this%water_table_old(icell) = this%water_table(icell)
  end subroutine sethead

  !> @brief Set infiltration
  !<
  subroutine setdatafinf(this, icell, finf)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: finf
    !
    if (this%landflag(icell) == 1) then
      this%finf_spec(icell) = finf
      this%finf(icell) = finf
    else
      this%finf_spec(icell) = DZERO
      this%finf(icell) = DZERO
    end if
    this%finf_rej(icell) = DZERO
    this%surf_infil(icell) = DZERO
    this%surf_infil_below(icell) = DZERO
  end subroutine setdatafinf

  !> @brief Set uzfarea using cellarea and areamult
  !<
  subroutine setdatauzfarea(this, icell, areamult)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: areamult
    !
    ! -- set uzf area
    this%uzfarea(icell) = this%cellarea(icell) * areamult
  end subroutine setdatauzfarea

  !> @brief Set unsaturated ET-related variables
  !<
  subroutine setdataet(this, icell, jbelow, pet, extdp)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: jbelow
    real(DP), intent(in) :: pet
    real(DP), intent(in) :: extdp
    ! -- local
    real(DP) :: thick
    !
    if (this%landflag(icell) == 1) then
      this%pet(icell) = pet
      this%gw_pet(icell) = pet
    else
      this%pet(icell) = DZERO
      this%gw_pet(icell) = DZERO
    end if
    thick = this%celtop(icell) - this%celbot(icell)
    this%ext_depth(icell) = extdp
    if (this%landflag(icell) > 0) then
      this%landtop(icell) = this%celtop(icell)
      this%pet_max(icell) = this%pet(icell)
    end if
    !
    ! -- set uz extinction depth
    if (this%landtop(icell) - this%ext_depth(icell) < this%celbot(icell)) then
      this%ext_depth_uz(icell) = thick
    else
      this%ext_depth_uz(icell) = this%celtop(icell) - &
                                 (this%landtop(icell) - this%ext_depth(icell))
    end if
    if (this%ext_depth_uz(icell) < DZERO) this%ext_depth_uz(icell) = DZERO
    if (this%ext_depth_uz(icell) > DEM7 .and. this%ext_depth(icell) < DEM7) &
      this%ext_depth(icell) = this%ext_depth_uz(icell)
    !
    ! -- set pet for underlying cell
    if (jbelow > 0) then
      this%landtop(jbelow) = this%landtop(icell)
      this%pet_max(jbelow) = this%pet_max(icell)
    end if
  end subroutine setdataet

  !> @brief Set extinction water content
  !<
  subroutine setdataetwc(this, icell, jbelow, extwc)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: jbelow
    real(DP), intent(in) :: extwc
    !
    ! -- set extinction water content
    this%theta_ext(icell) = extwc
    if (jbelow > 0) this%theta_ext(jbelow) = extwc
  end subroutine setdataetwc

  !> @brief Set variables for head-based unsaturated flow
  !<
  subroutine setdataetha(this, icell, jbelow, ha, hroot, rootact)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: jbelow
    real(DP), intent(in) :: ha
    real(DP), intent(in) :: hroot
    real(DP), intent(in) :: rootact
    !
    ! -- set variables
    this%air_entry(icell) = ha
    this%root_pot(icell) = hroot
    this%root_act(icell) = rootact
    if (jbelow > 0) then
      this%air_entry(jbelow) = ha
      this%root_pot(jbelow) = hroot
      this%root_act(jbelow) = rootact
    end if
  end subroutine setdataetha

  !> @brief Set variables to advance to new time step. nothing yet.
  !<
  subroutine advance(this, icell)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    !
    ! -- set variables
    this%surf_seep(icell) = DZERO
  end subroutine advance

  !> @brief Route water through one uzf cell and form its gwf equation terms
  !!
  !! With reset_state, the wave train is saved on entry and restored on exit, so
  !! a trial solution leaves no trace; only the final call of the time step, with
  !! reset_state false, commits the new waves through update_wav.
  !<
  subroutine solve(this, jbelow, icell, totfluxtot, ietflag, &
                   issflag, iseepflag, hgwf, qfrommvr, ierr, &
                   reset_state, trhs, thcof, deriv, watercontent)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: jbelow !< number of underlying uzf object or 0 if none
    integer(I4B), intent(in) :: icell !< number of this uzf object
    real(DP), intent(inout) :: totfluxtot !<
    integer(I4B), intent(in) :: ietflag !< et is off (0) or based one water content (1) or pressure (2)
    integer(I4B), intent(in) :: issflag !< steady state flag
    integer(I4B), intent(in) :: iseepflag !< discharge to land is active (1) or not (0)
    real(DP), intent(in) :: hgwf !< head for cell icell
    real(DP), intent(in) :: qfrommvr !< water inflow from mover
    integer(I4B), intent(inout) :: ierr !< flag indicating not enough waves
    logical, intent(in) :: reset_state !< flag indicating that waves should be reset after solution
    real(DP), intent(inout), optional :: trhs !< total uzf rhs contribution to GWF model
    real(DP), intent(inout), optional :: thcof !< total uzf hcof contribution to GWF model
    real(DP), intent(inout), optional :: deriv !< derivate term for contribution to GWF model
    real(DP), intent(inout), optional :: watercontent !< calculated water content
    ! -- local
    real(DP) :: test
    real(DP) :: scale
    real(DP) :: seep
    real(DP) :: finfact
    real(DP) :: derivfinf
    real(DP) :: trhsfinf
    real(DP) :: thcoffinf
    real(DP) :: trhsseep
    real(DP) :: thcofseep
    real(DP) :: deriv1
    real(DP) :: deriv2
    !
    ! -- initialize variables
    totfluxtot = DZERO
    trhsfinf = DZERO
    thcoffinf = DZERO
    trhsseep = DZERO
    thcofseep = DZERO
    this%finf_rej(icell) = DZERO
    this%surf_infil(icell) = this%finf(icell) + qfrommvr / this%uzfarea(icell)
    this%water_table(icell) = hgwf
    this%surf_seep(icell) = DZERO
    seep = DZERO
    finfact = DZERO
    this%et_uz(icell) = DZERO
    this%surf_infil_below(icell) = DZERO
    if (this%cell_below(icell) > 0) then
      this%finf(jbelow) = DZERO
    end if
    if (this%water_table(icell) < this%celbot(icell)) &
      this%water_table(icell) = this%celbot(icell)
    !
    ! -- initialize derivative variables
    deriv1 = DZERO
    deriv2 = DZERO
    derivfinf = DZERO
    !
    ! -- save wave states for resetting after iteration.
    if (reset_state) then
      call this%store_waves(icell, this%wavsav)
    end if
    !
    if (this%water_table(icell) > this%celtop(icell)) &
      this%water_table(icell) = this%celtop(icell)
    !
    ! -- add water from mover to applied infiltration.
    if (this%surf_infil(icell) > this%vks(icell)) then
      this%surf_infil(icell) = this%vks(icell)
    end if
    !
    ! -- saturation excess rejected infiltration
    if (this%landflag(icell) == 1) then
      call this%rejfinf(icell, deriv1, hgwf, trhsfinf, thcoffinf, finfact)
      this%surf_infil(icell) = finfact
    end if
    !
    ! -- calculate rejected infiltration
    this%finf_rej(icell) = this%finf(icell) + &
                         (qfrommvr / this%uzfarea(icell)) - this%surf_infil(icell)
    !
    ! -- calculate groundwater discharge
    if (iseepflag > 0 .and. this%landflag(icell) == 1) then
      call this%gwseep(icell, deriv2, scale, hgwf, trhsseep, thcofseep, seep)
      this%surf_seep(icell) = seep
    end if
    !
    ! -- route water through unsat zone, calc. storage change and recharge
    test = this%water_table(icell)
    if (this%water_table_old(icell) - test < -DEM15) test = this%water_table_old(icell)
    if (this%celtop(icell) - test > DEM15) then
      if (issflag == 0) then
        call this%routewaves(totfluxtot, delt, ietflag, icell, ierr)
        if (ierr > 0) then
          if (reset_state) &
            call this%load_waves(icell, this%wavsav)
          return
        end if
        call this%uz_rise(icell, totfluxtot)
        this%flux_to_wt(icell) = totfluxtot
        if (this%cell_below(icell) > 0) then
          call this%addrech(icell, jbelow, hgwf, trhsfinf, thcoffinf, &
                            derivfinf, delt)
        end if
      else
        this%flux_to_wt(icell) = this%surf_infil(icell) * delt
        totfluxtot = this%surf_infil(icell) * delt
      end if
      thcoffinf = DZERO
      trhsfinf = this%flux_to_wt(icell) * this%uzfarea(icell) / delt
      if (.not. reset_state) then
        call this%update_wav(icell, delt, issflag, 0)
      end if
    else
      this%flux_to_wt(icell) = this%surf_infil(icell) * delt
      totfluxtot = this%surf_infil(icell) * delt
      if (.not. reset_state) then
        call this%update_wav(icell, delt, issflag, 1)
      end if
    end if
    !
    ! -- If formulating, then these variables will be present
    if (present(deriv)) deriv = deriv1 + deriv2 + derivfinf
    if (present(trhs)) trhs = trhsfinf + trhsseep
    if (present(thcof)) thcof = thcoffinf + thcofseep
    !
    ! -- Assign water content prior to resetting waves
    if (present(watercontent)) then
      watercontent = this%get_wcnew(icell)
    end if
    !
    ! -- reset waves to previous state for next iteration
    if (reset_state) then
      call this%load_waves(icell, this%wavsav)
    end if
  end subroutine solve

  !> @brief Add recharge or infiltration to cells
  !<
  subroutine addrech(this, icell, jbelow, hgwf, trhs, thcof, deriv, delt)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: jbelow
    real(DP), intent(inout) :: trhs
    real(DP), intent(inout) :: thcof
    real(DP), intent(inout) :: deriv
    real(DP), intent(in) :: delt
    real(DP), intent(in) :: hgwf
    ! -- local
    real(DP) :: fcheck
    real(DP) :: x, scale, range
    !
    ! -- initialize
    range = DEM5
    deriv = DZERO
    thcof = DZERO
    trhs = this%uzfarea(icell) * this%flux_to_wt(icell) / delt
    if (this%flux_to_wt(icell) < DEM14) return
    scale = DONE
    !
    ! -- smoothly reduce flow between cells when head close to cell top
    x = hgwf - (this%celbot(icell) - range)
    call sSCurve(x, range, deriv, scale)
    deriv = this%uzfarea(icell) * deriv * this%flux_to_wt(icell) / delt
    this%finf(jbelow) = (DONE - scale) * this%flux_to_wt(icell) / delt
    fcheck = this%finf(jbelow) - this%vks(jbelow)
    !
    ! -- reduce flow between cells when vks is too small
    if (fcheck < DEM14) fcheck = DZERO
    this%finf(jbelow) = this%finf(jbelow) - fcheck
    this%surf_infil_below(icell) = this%finf(jbelow)
    this%flux_to_wt(icell) = scale * this%flux_to_wt(icell) + fcheck * delt
    trhs = this%uzfarea(icell) * this%flux_to_wt(icell) / delt
  end subroutine addrech

  !> @brief Reject applied infiltration due to low vks
  !<
  subroutine rejfinf(this, icell, deriv, hgwf, trhs, thcof, finfact)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(inout) :: deriv
    real(DP), intent(inout) :: finfact
    real(DP), intent(inout) :: thcof
    real(DP), intent(inout) :: trhs
    real(DP), intent(in) :: hgwf
    ! -- local
    real(DP) :: x, range, scale, q
    !
    range = this%surfdep(icell)
    q = this%surf_infil(icell)
    finfact = q
    trhs = finfact * this%uzfarea(icell)
    x = this%celtop(icell) - hgwf
    call sLinear(x, range, deriv, scale)
    deriv = -q * deriv * this%uzfarea(icell) * scale
    if (scale < DONE) then
      finfact = q * scale
      trhs = finfact * this%uzfarea(icell) * this%celtop(icell) / range
      thcof = finfact * this%uzfarea(icell) / range
    end if
  end subroutine rejfinf

  !> @brief Calculate groundwater discharge to land surface
  !<
  subroutine gwseep(this, icell, deriv, scale, hgwf, trhs, thcof, seep)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(inout) :: deriv
    real(DP), intent(inout) :: trhs
    real(DP), intent(inout) :: thcof
    real(DP), intent(inout) :: seep
    real(DP), intent(out) :: scale
    real(DP), intent(in) :: hgwf
    ! -- local
    real(DP) :: x, range, y, deriv1, d1, d2, Q
    !
    seep = DZERO
    deriv = DZERO
    deriv1 = DZERO
    d1 = DZERO
    d2 = DZERO
    scale = DZERO
    Q = this%uzfarea(icell) * this%vks(icell)
    range = this%surfdep(icell)
    x = hgwf - this%celtop(icell)
    call sCubicLinear(x, range, deriv1, y)
    scale = y
    seep = scale * Q * (hgwf - this%celtop(icell)) / range
    trhs = scale * Q * this%celtop(icell) / range
    thcof = -scale * Q / range
    d1 = -deriv1 * Q * x / range
    d2 = -scale * Q / range
    deriv = d1 + d2
    if (seep < DZERO) then
      seep = DZERO
      deriv = DZERO
      trhs = DZERO
      thcof = DZERO
    end if
  end subroutine gwseep

end module UzfCellGroupModule
