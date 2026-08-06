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

  !> @brief Subtract aet from pet to calculate residual et for gw
  !<
  subroutine setgwpet(this, icell)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    ! -- dummy
    real(DP) :: pet
    !
    pet = DZERO
    !
    ! -- reduce pet for gw by uzet
    pet = this%pet(icell) - this%et_uz(icell) / delt
    if (pet < DZERO) pet = DZERO
    this%gw_pet(icell) = pet
  end subroutine setgwpet

  !> @brief Subtract aet from pet to calculate residual et for deeper cells
  !<
  subroutine setbelowpet(this, icell, jbelow)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: jbelow
    ! -- dummy
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
  end subroutine setbelowpet

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

  !> @brief Calculate gwf et using residual uzf pet
  !<
  subroutine simgwet(this, igwetflag, icell, hgwf, trhs, thcof, det)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: igwetflag
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: hgwf
    real(DP), intent(inout) :: trhs
    real(DP), intent(inout) :: thcof
    real(DP), intent(inout) :: det
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
  end subroutine simgwet

  !> @brief Calculate recharge due to a rise in the gwf head
  !<
  subroutine uz_rise(this, icell, totfluxtot)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(inout) :: totfluxtot
    ! -- local
    real(DP) :: fm1, fm2, d1
    !
    ! -- additional recharge from a rising water table
    if (this%water_table(icell) - this%water_table_old(icell) > DEM30) then
      d1 = this%celtop(icell) - this%water_table_old(icell)
      fm1 = this%unsat_stor(icell, d1)
      d1 = this%celtop(icell) - this%water_table(icell)
      fm2 = this%unsat_stor(icell, d1)
      totfluxtot = totfluxtot + (fm1 - fm2)
    end if
  end subroutine uz_rise

  !> @brief Reset waves to default values at start of simulation
  !<
  subroutine setwaves(this, icell)
    ! -- dummy
    class(UzfCellGroupType) :: this
    ! -- local
    integer(I4B), intent(in) :: icell
    real(DP) :: bottom, top
    integer(I4B) :: jk
    real(DP) :: thick
    !
    ! -- initialize
    this%flux_to_wt(icell) = DZERO
    this%nwaves(icell) = 1
    this%wave_depth(:, icell) = DZERO
    thick = this%celtop(icell) - this%water_table(icell)
    do jk = 1, this%nwaves_max
      this%wave_theta(jk, icell) = this%theta_res(icell)
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
   this%wave_flux(1, icell) = this%vks(icell) * (top / bottom)**this%bc_eps(icell)
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
  end subroutine

  !> @brief Prepare and route waves over time step
  !<
  subroutine routewaves(this, totfluxtot, delt, ietflag, icell, ierr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    real(DP), intent(inout) :: totfluxtot
    real(DP), intent(in) :: delt
    integer(I4B), intent(in) :: ietflag
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(inout) :: ierr
    ! -- local
    real(DP) :: thick, thickold
    integer(I4B) :: idelt, iwav, ik
    !
    ! -- initialize
    this%flux_to_wt(icell) = DZERO
    this%et_uz(icell) = DZERO
    thick = this%celtop(icell) - this%water_table(icell)
    thickold = this%celtop(icell) - this%water_table_old(icell)
    !
    ! -- no uz, clear waves
    if (thickold < DZERO) then
      do iwav = 1, this%nwaves(icell)
        this%wave_theta(iwav, icell) = this%theta_res(icell)
        this%wave_depth(iwav, icell) = DZERO
        this%wave_speed(iwav, icell) = DZERO
        this%wave_flux(iwav, icell) = DZERO
      end do
      this%nwaves(icell) = 1
    end if
    idelt = 1
    do ik = 1, idelt
      call this%uzflow(thick, thickold, delt, ietflag, icell, ierr)
      if (ierr > 0) return
      totfluxtot = totfluxtot + this%flux_to_wt(icell)
    end do
  end subroutine routewaves

  !> @brief Move waves within a cell by shft positions
  !!
  !! Waves strt to stp are overwritten by the waves shft positions away.
  !! cntr is the loop increment; it must run the loop in the direction that
  !! keeps the source ahead of the destination.
  !<
  subroutine shift_waves(this, icell, shft, strt, stp, cntr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: shft
    integer(I4B), intent(in) :: strt
    integer(I4B), intent(in) :: stp
    integer(I4B), intent(in) :: cntr
    ! -- local
    integer(I4B) :: j
    !
    do j = strt, stp, cntr
      this%wave_theta(j, icell) = this%wave_theta(j + shft, icell)
      this%wave_depth(j, icell) = this%wave_depth(j + shft, icell)
      this%wave_flux(j, icell) = this%wave_flux(j + shft, icell)
      this%wave_speed(j, icell) = this%wave_speed(j + shft, icell)
    end do
  end subroutine shift_waves

  !> @brief Save the wave train of one cell into a snapshot
  !<
  subroutine store_waves(this, icell, store)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    type(UzfWaveStoreType), intent(inout) :: store
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
  end subroutine store_waves

  !> @brief Restore the wave train of one cell from a snapshot
  !<
  subroutine load_waves(this, icell, store)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    type(UzfWaveStoreType), intent(in) :: store
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
  end subroutine load_waves

  !> @brief Method of Characteristics solution for kinematic wave equation
  !<
  subroutine uzflow(this, thick, thickold, delt, ietflag, icell, ierr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    real(DP), intent(inout) :: thickold
    real(DP), intent(inout) :: thick
    real(DP), intent(in) :: delt
    integer(I4B), intent(in) :: ietflag
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(inout) :: ierr
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
          this%wave_speed(2, icell) = this%wave_flux(2, icell) / &
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
 dflux_surf = (this%surf_infil(icell) - this%wave_flux(this%nwaves(icell), icell))
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
  end subroutine uzflow

  !> @brief Calculate unit specific tolerances
  !<
  subroutine factors(eps_thick, eps_flux)
    ! -- dummy
    real(DP), intent(out) :: eps_thick
    real(DP), intent(out) :: eps_flux
    real(DP) :: factor1
    real(DP) :: factor2
    !
    ! calculate constants for uzflow
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
  end subroutine factors

  !> @brief Create and set trail waves
  !<
  subroutine trailwav(this, icell, ierr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(inout) :: ierr
    ! -- local
    real(DP) :: theta_surf, theta_step, ftrail, eps_m1
    real(DP) :: dtheta_inv
    real(DP) :: flux1, flux2, theta1, theta2
    real(DP) :: fnuminc
    integer(I4B) :: j, jj, jk, nwaves_m1
    !
    ! -- initialize
    eps_m1 = dble(this%bc_eps(icell)) - DONE
    dtheta_inv = DONE / (this%theta_sat(icell) - this%theta_res(icell))
    nwaves_m1 = this%nwaves(icell) - 1
    !
    ! -- initialize trailwaves
    theta_surf = (((this%surf_infil(icell) / this%vks(icell))** &
                   (DONE / this%bc_eps(icell))) * &
          (this%theta_sat(icell) - this%theta_res(icell))) + this%theta_res(icell)
    if (this%wave_theta(nwaves_m1, icell) - theta_surf > DEM9) then
      fnuminc = DZERO
      do jk = 1, this%ntrailwaves
        fnuminc = fnuminc + float(jk)
      end do
  theta_step = (this%wave_theta(nwaves_m1, icell) - theta_surf) / (fnuminc - DONE)
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
        this%wave_flux(j, icell) = &
       this%vks(icell) * (((this%wave_theta(j, icell) - this%theta_res(icell)) * &
                              dtheta_inv)**this%bc_eps(icell))
        theta2 = this%wave_theta(j - 1, icell)
        flux2 = this%wave_flux(j - 1, icell)
        flux1 = this%wave_flux(j, icell)
        theta1 = this%wave_theta(j, icell)
        this%wave_speed(j, icell) = leadspeed(theta1, theta2, flux1, flux2, &
                                   this%theta_sat(icell), this%theta_res(icell), &
                                              this%bc_eps(icell), this%vks(icell))
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
      this%wave_depth(this%nwaves(icell), icell) = DZERO
      this%wave_flux(this%nwaves(icell), icell) = &
        this%vks(icell) * (((this%wave_theta(this%nwaves(icell), icell) - &
                             this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
      this%wave_theta(this%nwaves(icell), icell) = theta_surf
      theta2 = this%wave_theta(this%nwaves(icell) - 1, icell)
      flux2 = this%wave_flux(this%nwaves(icell) - 1, icell)
      flux1 = this%wave_flux(this%nwaves(icell), icell)
      theta1 = this%wave_theta(this%nwaves(icell), icell)
      this%wave_speed(this%nwaves(icell), icell) = &
        leadspeed(theta1, theta2, flux1, flux2, this%theta_sat(icell), &
                  this%theta_res(icell), this%bc_eps(icell), this%vks(icell))
    end if
  end subroutine trailwav

  !> @brief Create a lead wave and route over time step
  !<
  subroutine leadwav(this, time, single_wave, trail_added, thetab, fluxb, &
                     dflux_surf, eps_flux, delt, icell)
    ! -- dummy
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
    ! -- local
    real(DP) :: dt_to_bottom, dt_first, fcheck
    real(DP) :: eps_m1, time_new, bottom, time_left
    real(DP) :: dtheta_inv, diff, flux_bot_prev
    real(DP) :: flux1, flux2, theta1, theta2, ftest
    integer(I4B) :: wave_exited, iremove, j, l
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
 (DONE / this%bc_eps(icell))) * (this%theta_sat(icell) - this%theta_res(icell))) &
          + this%theta_res(icell)
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
   dt_to_bottom = (this%wave_depth(1, icell) - this%wave_depth(2, icell)) / bottom
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
          l = j
          do while (j < this%nwaves(icell) + 1)
            if (this%wave_merges(j) == 1) then
              l = j
              theta2 = this%wave_theta(j, icell)
              flux2 = this%wave_flux(j, icell)
              if (j == 3) then
                flux1 = fluxb
                theta1 = thetab
              else
                flux1 = this%wave_flux(j - 2, icell)
                theta1 = this%wave_theta(j - 2, icell)
              end if
             this%wave_speed(j, icell) = leadspeed(theta1, theta2, flux1, flux2, &
                                                    this%theta_sat(icell), &
                                                    this%theta_res(icell), &
                                              this%bc_eps(icell), this%vks(icell))
              !
              ! -- update waves
              call this%shift_waves(icell, 1, l - 1, &
                                    this%nwaves(icell) - 1, 1)
              l = this%nwaves(icell) + 1
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
        this%flux_to_wt(icell) = this%flux_to_wt(icell) + flux_bot_prev * (time_new - time)
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
  end subroutine leadwav

  !> @brief Calculates waves speed from dflux/dtheta
  !<
  function leadspeed(theta1, theta2, flux1, flux2, thts, thtr, eps, vks)
    ! -- Return
    real(DP) :: leadspeed
    ! -- dummy
    real(DP), intent(in) :: theta1
    real(DP), intent(in) :: theta2
    real(DP), intent(in) :: flux1
    real(DP), intent(inout) :: flux2
    real(DP), intent(in) :: thts
    real(DP), intent(in) :: thtr
    real(DP), intent(in) :: eps
    real(DP), intent(in) :: vks
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
      fhold = DEM30
      if (comp3 > DEM30) fhold = (comp3 * thsrinv)**eps
      if (fhold < DEM30) fhold = DEM30
      leadspeed = epsfksths * (fhold**eps_m1)
    else
      leadspeed = (flux2 - flux1) / (theta2 - theta1)
    end if
    if (leadspeed < DEM30) leadspeed = DEM30
  end function leadspeed

  !> @brief Sums up mobile water over depth interval
  !<
  function unsat_stor(this, icell, d1)
    ! -- Return
    real(DP) :: unsat_stor
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(inout) :: d1
    ! -- local
    real(DP) :: fm
    integer(I4B) :: j, k, nwavm1, jj
    !
    fm = DZERO
    j = this%nwaves(icell) + 1
    k = this%nwaves(icell)
    nwavm1 = k - 1
    if (d1 > this%wave_depth(1, icell)) d1 = this%wave_depth(1, icell)
    !
    ! -- find deepest wave above depth d1, counter held as j
    do while (k > 0)
      if (this%wave_depth(k, icell) - d1 < -DEM30) j = k
      k = k - 1
    end do
    if (j > this%nwaves(icell)) then
      fm = fm + (this%wave_theta(this%nwaves(icell), icell) - this%theta_res(icell)) * d1
    elseif (this%nwaves(icell) > 1) then
      if (j > 1) then
        fm = fm + (this%wave_theta(j - 1, icell) - this%theta_res(icell)) &
             * (d1 - this%wave_depth(j, icell))
      end if
      do jj = j, nwavm1
        fm = fm + (this%wave_theta(jj, icell) - this%theta_res(icell)) &
             * (this%wave_depth(jj, icell) &
                - this%wave_depth(jj + 1, icell))
      end do
  fm = fm + (this%wave_theta(this%nwaves(icell), icell) - this%theta_res(icell)) &
           * (this%wave_depth(this%nwaves(icell), icell))
    else
      fm = fm + (this%wave_theta(1, icell) - this%theta_res(icell)) * d1
    end if
    unsat_stor = fm
  end function unsat_stor

  !> @brief Update to new state of uz at end of time step
  !<
  subroutine update_wav(this, icell, delt, iss, itest)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: itest
    integer(I4B), intent(in) :: iss
    real(DP), intent(in) :: delt
    ! -- local
    real(DP) :: bot, depthsave, top
    real(DP) :: thick, dtheta_inv
    integer(I4B) :: nwavhld, k, j
    !
    bot = this%water_table(icell)
    top = this%celtop(icell)
    thick = top - bot
    nwavhld = this%nwaves(icell)
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
        j = 0
        k = this%nwaves(icell)
        do while (k > 0)
          if (this%wave_depth(k, icell) - thick < -DEM30) j = k
          k = k - 1
        end do
        this%wave_depth(1, icell) = thick
        if (j > 1) then
          this%wave_speed(1, icell) = DZERO
          this%nwaves(icell) = this%nwaves(icell) - j + 2
          this%wave_theta(1, icell) = this%wave_theta(j - 1, icell)
          this%wave_flux(1, icell) = this%wave_flux(j - 1, icell)
          if (j > 2) call this%shift_waves(icell, j - 2, 2, &
                                           nwavhld - (j - 2), 1)
        elseif (j == 0) then
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
  end subroutine update_wav

  !> @brief Remove water from uz due to et
  !<
  subroutine uzet(this, icell, delt, ietflag, ierr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: delt
    integer(I4B), intent(in) :: ietflag
    integer(I4B), intent(inout) :: ierr
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
    real(DP) :: petsub
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: jhold
    integer(I4B) :: jk
    integer(I4B) :: kj
    integer(I4B) :: kk
    integer(I4B) :: numadd
    integer(I4B) :: k
    integer(I4B) :: nwv
    integer(I4B) :: itest
    !
    ! -- initialize
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
 if ((this%wave_theta(1, icell) - thetaout) > this%theta_res(icell) + extwc1) then
          this%wave_theta(1, icell) = this%wave_theta(1, icell) - thetaout
          this%wave_flux(1, icell) = &
            this%vks(icell) * (((this%wave_theta(1, icell) - &
                         this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
        else if (this%wave_theta(1, icell) > this%theta_res(icell) + extwc1) then
          this%wave_theta(1, icell) = this%theta_res(icell) + extwc1
          this%wave_flux(1, icell) = &
            this%vks(icell) * (((this%wave_theta(1, icell) - &
                         this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
        end if
        !
        ! -- all waves shallower than extinction depth
      else if (this%nwaves(icell) > 1 .AND. &
       this%wave_depth(this%nwaves(icell), icell) > this%ext_depth_uz(icell)) then
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
            this%theta_res(icell) + extwc1) then
          this%wave_theta(this%nwaves(icell) + 1, icell) = &
            this%wave_theta(this%nwaves(icell), icell) - thetaout
          numadd = 1
        else if (this%wave_theta(this%nwaves(icell), icell) > &
                 this%theta_res(icell) + extwc1) then
   this%wave_theta(this%nwaves(icell) + 1, icell) = this%theta_res(icell) + extwc1
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
            ! -- too many waves error, deallocate temp arrays and return
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
 if ((this%wave_theta(1, icell) - thetaout) > this%theta_res(icell) + extwc1) then
          if (thetaout > DEM30) then
            this%wave_theta(2, icell) = this%wave_theta(1, icell) - thetaout
            this%wave_flux(2, icell) = &
       this%vks(icell) * (((this%wave_theta(2, icell) - this%theta_res(icell)) * &
                                  dtheta_inv)**this%bc_eps(icell))
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
        else if (this%wave_theta(1, icell) > this%theta_res(icell) + extwc1) then
          if (thetaout > DEM30) then
            this%wave_theta(2, icell) = this%theta_res(icell) + extwc1
            this%wave_flux(2, icell) = &
              this%vks(icell) * (((this%wave_theta(2, icell) - &
                         this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
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
          j = 2
          jk = 0
          !
          ! -- locate extinction depth between waves
          do while (jk == 0)
            diff = this%wave_depth(j, icell) - this%ext_depth_uz(icell)
            if (diff > dzero) then
              j = j + 1
            else
              jk = 1
            end if
          end do
          kk = j
          if (this%wave_theta(j, icell) > this%theta_res(icell) + extwc1) then
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
                                    this%nwaves(icell) + 1, j, -1)
              this%wave_depth(j, icell) = this%ext_depth_uz(icell)
              this%nwaves(icell) = this%nwaves(icell) + 1
              if (this%nwaves(icell) > this%nwaves_max) then
                !
                ! -- too many waves error
                ierr = 1
                goto 500
              end if
            end if
            kk = j
          else
            jhold = this%nwaves(icell)
            i = j + 1
            do while (i < this%nwaves(icell))
              if (this%wave_theta(i, icell) > this%theta_res(icell) + extwc1) then
                jhold = i
                i = this%nwaves(icell) + 1
              end if
              i = i + 1
            end do
            j = jhold
            kk = jhold
          end if
        else
          kk = 1
        end if
        !
        ! -- all waves above extinction depth
        do while (kk <= this%nwaves(icell))
          if (ietflag == 2) then
            tho = this%wave_theta(kk, icell)
            fktho = this%wave_flux(kk, icell)
            hcap = this%caph(icell, tho)
            thetaout = this%rate_et_z(icell, factor, fktho, hcap)
          end if
          if (this%wave_theta(kk, icell) > this%theta_res(icell) + extwc1) then
            if (this%wave_theta(kk, icell) - thetaout > &
                this%theta_res(icell) + extwc1) then
              this%wave_theta(kk, icell) = this%wave_theta(kk, icell) - thetaout
        else if (this%wave_theta(kk, icell) > this%theta_res(icell) + extwc1) then
              this%wave_theta(kk, icell) = this%theta_res(icell) + extwc1
            end if
            if (kk == 1) then
              this%wave_flux(kk, icell) = &
                this%vks(icell) * &
                (((this%wave_theta(kk, icell) - &
                   this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell))
            end if
            if (kk > 1) then
              flux1 = &
                this%vks(icell) * ((this%wave_theta(kk - 1, icell) - &
                          this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell)
              flux2 = &
                this%vks(icell) * ((this%wave_theta(kk, icell) - &
                          this%theta_res(icell)) * dtheta_inv)**this%bc_eps(icell)
              this%wave_flux(kk, icell) = flux2
              theta2 = this%wave_theta(kk, icell)
              theta1 = this%wave_theta(kk - 1, icell)
            this%wave_speed(kk, icell) = leadspeed(theta1, theta2, flux1, flux2, &
                                                     this%theta_sat(icell), &
                                                     this%theta_res(icell), &
                                              this%bc_eps(icell), this%vks(icell))
            end if
          end if
          kk = kk + 1
        end do
      end if
      !
      ! -- calculate aet
      kj = 1
      do while (kj <= this%nwaves(icell) - 1)
 if (abs(this%wave_theta(kj, icell) - this%wave_theta(kj + 1, icell)) < DEM6) then
          call this%shift_waves(icell, 1, kj + 1, &
                                this%nwaves(icell) - 1, 1)
          kj = kj - 1
          this%nwaves(icell) = this%nwaves(icell) - 1
        end if
        kj = kj + 1
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
  end subroutine uzet

  !> @brief Calculate capillary pressure head from B-C equation
  !<
  function caph(this, icell, tho)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: tho
    ! -- local
    real(DP) :: caph, lambda, star
    !
    caph = -DEM6
    star = (tho - this%theta_res(icell)) / (this%theta_sat(icell) - this%theta_res(icell))
    if (star < DEM15) star = DEM15
    lambda = DTWO / (this%bc_eps(icell) - DTHREE)
    if (star > DEM15) then
      if (tho - this%theta_sat(icell) < DEM15) then
        caph = this%air_entry(icell) * star**(-DONE / lambda)
      else
        caph = DZERO
      end if
    end if
  end function caph

  !> @brief Calculate capillary pressure-based uz et
  !<
  function rate_et_z(this, icell, factor, fktho, h)
    ! -- Return
    real(DP) :: rate_et_z
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: factor, fktho, h
    !
    rate_et_z = factor * fktho * (h - this%root_pot(icell))
    if (rate_et_z < DZERO) rate_et_z = DZERO
  end function rate_et_z

  !> @brief Determine the water content at a specific depth
  !!
  !! Because UZF-calculated waves are internal to UZF objects, different water
  !! contents exists at different depths.
  !<
  function get_water_content_at_depth(this, icell, depth) result(theta_at_depth)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell !< uzf cell containing depth
    real(DP), intent(in) :: depth !< depth within the cell
    ! -- return
    real(DP) :: theta_at_depth
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
  end function get_water_content_at_depth

  !> @brief Calculate and return the cell-based water content value
  !<
  function get_wcnew(this, icell) result(watercontent)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell !< uzf cell containing depth
    ! -- return
    real(DP) :: watercontent
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
  end function get_wcnew

end module UzfCellGroupModule
