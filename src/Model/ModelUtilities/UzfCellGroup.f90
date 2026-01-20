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

  type :: UzfCellGroupType

    integer(I4B) :: imem_manager
    real(DP), pointer, dimension(:), contiguous :: thtr => null()
    real(DP), pointer, dimension(:), contiguous :: thts => null()
    real(DP), pointer, dimension(:), contiguous :: thti => null()
    real(DP), pointer, dimension(:), contiguous :: eps => null()
    real(DP), pointer, dimension(:), contiguous :: extwc => null()
    real(DP), pointer, dimension(:), contiguous :: ha => null()
    real(DP), pointer, dimension(:), contiguous :: hroot => null()
    real(DP), pointer, dimension(:), contiguous :: rootact => null()
    real(DP), pointer, dimension(:), contiguous :: etact => null()
    real(DP), dimension(:, :), pointer, contiguous :: uzspst => null()
    real(DP), dimension(:, :), pointer, contiguous :: uzthst => null()
    real(DP), dimension(:, :), pointer, contiguous :: uzflst => null()
    real(DP), dimension(:, :), pointer, contiguous :: uzdpst => null()
    integer(I4B), pointer, dimension(:), contiguous :: nwavst => null()
    real(DP), pointer, dimension(:), contiguous :: totflux => null()
    integer(I4B), pointer, dimension(:), contiguous :: nwav => null()
    integer(I4B), pointer, dimension(:), contiguous :: ntrail => null()
    real(DP), pointer, dimension(:), contiguous :: sinf => null()
    real(DP), pointer, dimension(:), contiguous :: finf => null()
    real(DP), pointer, dimension(:), contiguous :: pet => null()
    real(DP), pointer, dimension(:), contiguous :: petmax => null()
    real(DP), pointer, dimension(:), contiguous :: extdp => null()
    real(DP), pointer, dimension(:), contiguous :: extdpuz => null()
    real(DP), pointer, dimension(:), contiguous :: finf_rej => null()
    real(DP), pointer, dimension(:), contiguous :: gwet => null()
    real(DP), pointer, dimension(:), contiguous :: uzfarea => null()
    real(DP), pointer, dimension(:), contiguous :: cellarea => null()
    real(DP), pointer, dimension(:), contiguous :: celtop => null()
    real(DP), pointer, dimension(:), contiguous :: celbot => null()
    real(DP), pointer, dimension(:), contiguous :: landtop => null()
    real(DP), pointer, dimension(:), contiguous :: watab => null()
    real(DP), pointer, dimension(:), contiguous :: watabold => null()
    real(DP), pointer, dimension(:), contiguous :: vks => null()
    real(DP), pointer, dimension(:), contiguous :: surfdep => null()
    real(DP), pointer, dimension(:), contiguous :: surflux => null()
    real(DP), pointer, dimension(:), contiguous :: surfluxbelow => null()
    real(DP), pointer, dimension(:), contiguous :: surfseep => null()
    real(DP), pointer, dimension(:), contiguous :: gwpet => null()
    integer(I4B), pointer, dimension(:), contiguous :: landflag => null()
    integer(I4B), pointer, dimension(:), contiguous :: ivertcon => null()

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
    procedure :: wave_shift
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
    integer(I4B), intent(in) :: nwav
    integer(I4B), intent(in) :: ncells
    character(len=*), intent(in), optional :: memory_path
    ! -- local
    integer(I4B) :: icell
    !
    ! -- Use mem_allocate if memory path is passed in, otherwise it's a temp object
    if (present(memory_path)) then
      this%imem_manager = 1
      call mem_allocate(this%uzdpst, nwav, ncells, 'UZDPST', memory_path)
      call mem_allocate(this%uzthst, nwav, ncells, 'UZTHST', memory_path)
      call mem_allocate(this%uzflst, nwav, ncells, 'UZFLST', memory_path)
      call mem_allocate(this%uzspst, nwav, ncells, 'UZSPST', memory_path)
      call mem_allocate(this%nwavst, ncells, 'NWAVST', memory_path)
      call mem_allocate(this%thtr, ncells, 'THTR', memory_path)
      call mem_allocate(this%thts, ncells, 'THTS', memory_path)
      call mem_allocate(this%thti, ncells, 'THTI', memory_path)
      call mem_allocate(this%eps, ncells, 'EPS', memory_path)
      call mem_allocate(this%ha, ncells, 'HA', memory_path)
      call mem_allocate(this%hroot, ncells, 'HROOT', memory_path)
      call mem_allocate(this%rootact, ncells, 'ROOTACT', memory_path)
      call mem_allocate(this%extwc, ncells, 'EXTWC', memory_path)
      call mem_allocate(this%etact, ncells, 'ETACT', memory_path)
      call mem_allocate(this%nwav, ncells, 'NWAV', memory_path)
      call mem_allocate(this%ntrail, ncells, 'NTRAIL', memory_path)
      call mem_allocate(this%totflux, ncells, 'TOTFLUX', memory_path)
      call mem_allocate(this%sinf, ncells, 'SINF', memory_path)
      call mem_allocate(this%finf, ncells, 'FINF', memory_path)
      call mem_allocate(this%finf_rej, ncells, 'FINF_REJ', memory_path)
      call mem_allocate(this%gwet, ncells, 'GWET', memory_path)
      call mem_allocate(this%uzfarea, ncells, 'UZFAREA', memory_path)
      call mem_allocate(this%cellarea, ncells, 'CELLAREA', memory_path)
      call mem_allocate(this%celtop, ncells, 'CELTOP', memory_path)
      call mem_allocate(this%celbot, ncells, 'CELBOT', memory_path)
      call mem_allocate(this%landtop, ncells, 'LANDTOP', memory_path)
      call mem_allocate(this%watab, ncells, 'WATAB', memory_path)
      call mem_allocate(this%watabold, ncells, 'WATABOLD', memory_path)
      call mem_allocate(this%surfdep, ncells, 'SURFDEP', memory_path)
      call mem_allocate(this%vks, ncells, 'VKS', memory_path)
      call mem_allocate(this%surflux, ncells, 'SURFLUX', memory_path)
      call mem_allocate(this%surfluxbelow, ncells, 'SURFLUXBELOW', memory_path)
      call mem_allocate(this%surfseep, ncells, 'SURFSEEP', memory_path)
      call mem_allocate(this%gwpet, ncells, 'GWPET', memory_path)
      call mem_allocate(this%pet, ncells, 'PET', memory_path)
      call mem_allocate(this%petmax, ncells, 'PETMAX', memory_path)
      call mem_allocate(this%extdp, ncells, 'EXTDP', memory_path)
      call mem_allocate(this%extdpuz, ncells, 'EXTDPUZ', memory_path)
      call mem_allocate(this%landflag, ncells, 'LANDFLAG', memory_path)
      call mem_allocate(this%ivertcon, ncells, 'IVERTCON', memory_path)
    else
      this%imem_manager = 0
      allocate (this%uzdpst(nwav, ncells))
      allocate (this%uzthst(nwav, ncells))
      allocate (this%uzflst(nwav, ncells))
      allocate (this%uzspst(nwav, ncells))
      allocate (this%nwavst(ncells))
      allocate (this%thtr(ncells))
      allocate (this%thts(ncells))
      allocate (this%thti(ncells))
      allocate (this%eps(ncells))
      allocate (this%ha(ncells))
      allocate (this%hroot(ncells))
      allocate (this%rootact(ncells))
      allocate (this%extwc(ncells))
      allocate (this%etact(ncells))
      allocate (this%nwav(ncells))
      allocate (this%ntrail(ncells))
      allocate (this%totflux(ncells))
      allocate (this%sinf(ncells))
      allocate (this%finf(ncells))
      allocate (this%finf_rej(ncells))
      allocate (this%gwet(ncells))
      allocate (this%uzfarea(ncells))
      allocate (this%cellarea(ncells))
      allocate (this%celtop(ncells))
      allocate (this%celbot(ncells))
      allocate (this%landtop(ncells))
      allocate (this%watab(ncells))
      allocate (this%watabold(ncells))
      allocate (this%surfdep(ncells))
      allocate (this%vks(ncells))
      allocate (this%surflux(ncells))
      allocate (this%surfluxbelow(ncells))
      allocate (this%surfseep(ncells))
      allocate (this%gwpet(ncells))
      allocate (this%pet(ncells))
      allocate (this%petmax(ncells))
      allocate (this%extdp(ncells))
      allocate (this%extdpuz(ncells))
      allocate (this%landflag(ncells))
      allocate (this%ivertcon(ncells))
    end if
    do icell = 1, ncells
      this%uzdpst(:, icell) = DZERO
      this%uzthst(:, icell) = DZERO
      this%uzflst(:, icell) = DZERO
      this%uzspst(:, icell) = DZERO
      this%nwavst(icell) = 1
      this%thtr(icell) = DZERO
      this%thts(icell) = DZERO
      this%thti(icell) = DZERO
      this%eps(icell) = DZERO
      this%ha(icell) = DZERO
      this%hroot(icell) = DZERO
      this%rootact(icell) = DZERO
      this%extwc(icell) = DZERO
      this%etact(icell) = DZERO
      this%nwav(icell) = nwav
      this%ntrail(icell) = 0
      this%totflux(icell) = DZERO
      this%sinf(icell) = DZERO
      this%finf(icell) = DZERO
      this%finf_rej(icell) = DZERO
      this%gwet(icell) = DZERO
      this%uzfarea(icell) = DZERO
      this%cellarea(icell) = DZERO
      this%celtop(icell) = DZERO
      this%celbot(icell) = DZERO
      this%landtop(icell) = DZERO
      this%watab(icell) = DZERO
      this%watabold(icell) = DZERO
      this%surfdep(icell) = DZERO
      this%vks(icell) = DZERO
      this%surflux(icell) = DZERO
      this%surfluxbelow(icell) = DZERO
      this%surfseep(icell) = DZERO
      this%gwpet(icell) = DZERO
      this%pet(icell) = DZERO
      this%petmax(icell) = DZERO
      this%extdp(icell) = DZERO
      this%extdpuz(icell) = DZERO
      this%landflag(icell) = 0
      this%ivertcon(icell) = 0
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
    ! -- deallocate based on whether or not memory manager was used
    if (this%imem_manager == 0) then
      deallocate (this%uzdpst)
      deallocate (this%uzthst)
      deallocate (this%uzflst)
      deallocate (this%uzspst)
      deallocate (this%nwavst)
      deallocate (this%thtr)
      deallocate (this%thts)
      deallocate (this%thti)
      deallocate (this%eps)
      deallocate (this%ha)
      deallocate (this%hroot)
      deallocate (this%rootact)
      deallocate (this%extwc)
      deallocate (this%etact)
      deallocate (this%nwav)
      deallocate (this%ntrail)
      deallocate (this%totflux)
      deallocate (this%sinf)
      deallocate (this%finf)
      deallocate (this%finf_rej)
      deallocate (this%gwet)
      deallocate (this%uzfarea)
      deallocate (this%cellarea)
      deallocate (this%celtop)
      deallocate (this%celbot)
      deallocate (this%landtop)
      deallocate (this%watab)
      deallocate (this%watabold)
      deallocate (this%surfdep)
      deallocate (this%vks)
      deallocate (this%surflux)
      deallocate (this%surfluxbelow)
      deallocate (this%surfseep)
      deallocate (this%gwpet)
      deallocate (this%pet)
      deallocate (this%petmax)
      deallocate (this%extdp)
      deallocate (this%extdpuz)
      deallocate (this%landflag)
      deallocate (this%ivertcon)
    else
      call mem_deallocate(this%uzdpst)
      call mem_deallocate(this%uzthst)
      call mem_deallocate(this%uzflst)
      call mem_deallocate(this%uzspst)
      call mem_deallocate(this%nwavst)
      call mem_deallocate(this%thtr)
      call mem_deallocate(this%thts)
      call mem_deallocate(this%thti)
      call mem_deallocate(this%eps)
      call mem_deallocate(this%ha)
      call mem_deallocate(this%hroot)
      call mem_deallocate(this%rootact)
      call mem_deallocate(this%extwc)
      call mem_deallocate(this%etact)
      call mem_deallocate(this%nwav)
      call mem_deallocate(this%ntrail)
      call mem_deallocate(this%totflux)
      call mem_deallocate(this%sinf)
      call mem_deallocate(this%finf)
      call mem_deallocate(this%finf_rej)
      call mem_deallocate(this%gwet)
      call mem_deallocate(this%uzfarea)
      call mem_deallocate(this%cellarea)
      call mem_deallocate(this%celtop)
      call mem_deallocate(this%celbot)
      call mem_deallocate(this%landtop)
      call mem_deallocate(this%watab)
      call mem_deallocate(this%watabold)
      call mem_deallocate(this%surfdep)
      call mem_deallocate(this%vks)
      call mem_deallocate(this%surflux)
      call mem_deallocate(this%surfluxbelow)
      call mem_deallocate(this%surfseep)
      call mem_deallocate(this%gwpet)
      call mem_deallocate(this%pet)
      call mem_deallocate(this%petmax)
      call mem_deallocate(this%extdp)
      call mem_deallocate(this%extdpuz)
      call mem_deallocate(this%landflag)
      call mem_deallocate(this%ivertcon)
    end if
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
    this%ivertcon(icell) = ivertcon
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
    this%thtr(icell) = thtr
    this%thts(icell) = thts
    this%thti(icell) = thti
    this%eps(icell) = eps
    this%ntrail(icell) = ntrail
    this%pet(icell) = DZERO
    this%extdp(icell) = DZERO
    this%extwc(icell) = DZERO
    this%ha(icell) = DZERO
    this%hroot(icell) = DZERO
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
    this%watab(icell) = this%celbot(icell)
    if (hgwf > this%celbot(icell)) this%watab(icell) = hgwf
    if (this%watab(icell) > this%celtop(icell)) &
      this%watab(icell) = this%celtop(icell)
    this%watabold(icell) = this%watab(icell)
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
      this%sinf(icell) = finf
      this%finf(icell) = finf
    else
      this%sinf(icell) = DZERO
      this%finf(icell) = DZERO
    end if
    this%finf_rej(icell) = DZERO
    this%surflux(icell) = DZERO
    this%surfluxbelow(icell) = DZERO
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
      this%gwpet(icell) = pet
    else
      this%pet(icell) = DZERO
      this%gwpet(icell) = DZERO
    end if
    thick = this%celtop(icell) - this%celbot(icell)
    this%extdp(icell) = extdp
    if (this%landflag(icell) > 0) then
      this%landtop(icell) = this%celtop(icell)
      this%petmax(icell) = this%pet(icell)
    end if
    !
    ! -- set uz extinction depth
    if (this%landtop(icell) - this%extdp(icell) < this%celbot(icell)) then
      this%extdpuz(icell) = thick
    else
      this%extdpuz(icell) = this%celtop(icell) - &
                            (this%landtop(icell) - this%extdp(icell))
    end if
    if (this%extdpuz(icell) < DZERO) this%extdpuz(icell) = DZERO
    if (this%extdpuz(icell) > DEM7 .and. this%extdp(icell) < DEM7) &
      this%extdp(icell) = this%extdpuz(icell)
    !
    ! -- set pet for underlying cell
    if (jbelow > 0) then
      this%landtop(jbelow) = this%landtop(icell)
      this%petmax(jbelow) = this%petmax(icell)
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
    pet = this%pet(icell) - this%etact(icell) / delt
    if (pet < DZERO) pet = DZERO
    this%gwpet(icell) = pet
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
    if (this%extdpuz(jbelow) > DEM3) then
      pet = this%pet(icell) - this%etact(icell) / delt - &
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
    this%extwc(icell) = extwc
    if (jbelow > 0) this%extwc(jbelow) = extwc
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
    this%ha(icell) = ha
    this%hroot(icell) = hroot
    this%rootact(icell) = rootact
    if (jbelow > 0) then
      this%ha(jbelow) = ha
      this%hroot(jbelow) = hroot
      this%rootact(jbelow) = rootact
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
    this%surfseep(icell) = DZERO
  end subroutine advance

  !> @brief Formulate the unsaturated flow object, calculate terms for gwf
  !! equation
  !<
  subroutine solve(this, thiswork, jbelow, icell, totfluxtot, ietflag, &
                   issflag, iseepflag, hgwf, qfrommvr, ierr, &
                   reset_state, trhs, thcof, deriv, watercontent)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(UzfCellGroupType) :: this
    type(UzfCellGroupType) :: thiswork !< work object for resetting wave state
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
    this%surflux(icell) = this%finf(icell) + qfrommvr / this%uzfarea(icell)
    this%watab(icell) = hgwf
    this%surfseep(icell) = DZERO
    seep = DZERO
    finfact = DZERO
    this%etact(icell) = DZERO
    this%surfluxbelow(icell) = DZERO
    if (this%ivertcon(icell) > 0) then
      this%finf(jbelow) = DZERO
    end if
    if (this%watab(icell) < this%celbot(icell)) &
      this%watab(icell) = this%celbot(icell)
    !
    ! -- initialize derivative variables
    deriv1 = DZERO
    deriv2 = DZERO
    derivfinf = DZERO
    !
    ! -- save wave states for resetting after iteration.
    if (reset_state) then
      call thiswork%wave_shift(this, 1, icell, 0, 1, this%nwavst(icell), 1)
    end if
    !
    if (this%watab(icell) > this%celtop(icell)) &
      this%watab(icell) = this%celtop(icell)
    !
    ! -- add water from mover to applied infiltration.
    if (this%surflux(icell) > this%vks(icell)) then
      this%surflux(icell) = this%vks(icell)
    end if
    !
    ! -- saturation excess rejected infiltration
    if (this%landflag(icell) == 1) then
      call this%rejfinf(icell, deriv1, hgwf, trhsfinf, thcoffinf, finfact)
      this%surflux(icell) = finfact
    end if
    !
    ! -- calculate rejected infiltration
    this%finf_rej(icell) = this%finf(icell) + &
                           (qfrommvr / this%uzfarea(icell)) - this%surflux(icell)
    !
    ! -- calculate groundwater discharge
    if (iseepflag > 0 .and. this%landflag(icell) == 1) then
      call this%gwseep(icell, deriv2, scale, hgwf, trhsseep, thcofseep, seep)
      this%surfseep(icell) = seep
    end if
    !
    ! -- route water through unsat zone, calc. storage change and recharge
    test = this%watab(icell)
    if (this%watabold(icell) - test < -DEM15) test = this%watabold(icell)
    if (this%celtop(icell) - test > DEM15) then
      if (issflag == 0) then
        call this%routewaves(totfluxtot, delt, ietflag, icell, ierr)
        if (ierr > 0) return
        call this%uz_rise(icell, totfluxtot)
        this%totflux(icell) = totfluxtot
        if (this%ivertcon(icell) > 0) then
          call this%addrech(icell, jbelow, hgwf, trhsfinf, thcoffinf, &
                            derivfinf, delt)
        end if
      else
        this%totflux(icell) = this%surflux(icell) * delt
        totfluxtot = this%surflux(icell) * delt
      end if
      thcoffinf = DZERO
      trhsfinf = this%totflux(icell) * this%uzfarea(icell) / delt
      if (.not. reset_state) then
        call this%update_wav(icell, delt, issflag, 0)
      end if
    else
      this%totflux(icell) = this%surflux(icell) * delt
      totfluxtot = this%surflux(icell) * delt
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
      call this%wave_shift(thiswork, icell, 1, 0, 1, thiswork%nwavst(1), 1)
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
    trhs = this%uzfarea(icell) * this%totflux(icell) / delt
    if (this%totflux(icell) < DEM14) return
    scale = DONE
    !
    ! -- smoothly reduce flow between cells when head close to cell top
    x = hgwf - (this%celbot(icell) - range)
    call sSCurve(x, range, deriv, scale)
    deriv = this%uzfarea(icell) * deriv * this%totflux(icell) / delt
    this%finf(jbelow) = (DONE - scale) * this%totflux(icell) / delt
    fcheck = this%finf(jbelow) - this%vks(jbelow)
    !
    ! -- reduce flow between cells when vks is too small
    if (fcheck < DEM14) fcheck = DZERO
    this%finf(jbelow) = this%finf(jbelow) - fcheck
    this%surfluxbelow(icell) = this%finf(jbelow)
    this%totflux(icell) = scale * this%totflux(icell) + fcheck * delt
    trhs = this%uzfarea(icell) * this%totflux(icell) / delt
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
    q = this%surflux(icell)
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
    x = this%extdp(icell)
    c = this%gwpet(icell)
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
    if (this%watab(icell) - this%watabold(icell) > DEM30) then
      d1 = this%celtop(icell) - this%watabold(icell)
      fm1 = this%unsat_stor(icell, d1)
      d1 = this%celtop(icell) - this%watab(icell)
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
    this%totflux(icell) = DZERO
    this%nwavst(icell) = 1
    this%uzdpst(:, icell) = DZERO
    thick = this%celtop(icell) - this%watab(icell)
    do jk = 1, this%nwav(icell)
      this%uzthst(jk, icell) = this%thtr(icell)
    end do
    !
    ! -- initialize waves for first stress period
    if (thick > DZERO) then
      this%uzdpst(1, icell) = thick
      this%uzthst(1, icell) = this%thti(icell)
      top = this%uzthst(1, icell) - this%thtr(icell)
      if (top < DZERO) top = DZERO
      bottom = this%thts(icell) - this%thtr(icell)
      if (bottom < DZERO) bottom = DZERO
      this%uzflst(1, icell) = this%vks(icell) * (top / bottom)**this%eps(icell)
      if (this%uzthst(1, icell) < this%thtr(icell)) &
        this%uzthst(1, icell) = this%thtr(icell)
      !
      ! -- calculate water stored in the unsaturated zone
      if (top > DZERO) then
        this%uzspst(1, icell) = DZERO
      else
        this%uzflst(1, icell) = DZERO
        this%uzspst(1, icell) = DZERO
      end if
      !
      ! no unsaturated zone
    else
      this%uzflst(1, icell) = DZERO
      this%uzdpst(1, icell) = DZERO
      this%uzspst(1, icell) = DZERO
      this%uzthst(1, icell) = this%thtr(icell)
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
    this%totflux(icell) = DZERO
    this%etact(icell) = DZERO
    thick = this%celtop(icell) - this%watab(icell)
    thickold = this%celtop(icell) - this%watabold(icell)
    !
    ! -- no uz, clear waves
    if (thickold < DZERO) then
      do iwav = 1, 6
        this%uzthst(iwav, icell) = this%thtr(icell)
        this%uzdpst(iwav, icell) = DZERO
        this%uzspst(iwav, icell) = DZERO
        this%uzflst(iwav, icell) = DZERO
        this%nwavst(icell) = 1
      end do
    end if
    idelt = 1
    do ik = 1, idelt
      call this%uzflow(thick, thickold, delt, ietflag, icell, ierr)
      if (ierr > 0) return
      totfluxtot = totfluxtot + this%totflux(icell)
    end do
  end subroutine routewaves

  !> @brief Copy waves or shift waves in arrays
  !<
  subroutine wave_shift(this, this2, icell, icell2, shft, strt, stp, cntr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    type(UzfCellGroupType) :: this2
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(in) :: icell2
    integer(I4B), intent(in) :: shft
    integer(I4B), intent(in) :: strt
    integer(I4B), intent(in) :: stp
    integer(I4B), intent(in) :: cntr
    ! -- local
    integer(I4B) :: j
    !
    ! -- copy waves from one uzf cell group to another
    do j = strt, stp, cntr
      this%uzthst(j, icell) = this2%uzthst(j + shft, icell2)
      this%uzdpst(j, icell) = this2%uzdpst(j + shft, icell2)
      this%uzflst(j, icell) = this2%uzflst(j + shft, icell2)
      this%uzspst(j, icell) = this2%uzspst(j + shft, icell2)
    end do
    this%nwavst(icell) = this2%nwavst(icell2)
  end subroutine

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
    real(DP) :: ffcheck, time, feps1, feps2
    real(DP) :: thetadif, thetab, fluxb, oldsflx
    integer(I4B) :: itrailflg, itester
    !
    time = DZERO
    this%totflux(icell) = DZERO
    itrailflg = 0
    oldsflx = this%uzflst(this%nwavst(icell), icell)
    call factors(feps1, feps2)
    !
    ! -- check for falling or rising water table
    if ((thick - thickold) > feps1) then
      thetadif = abs(this%uzthst(1, icell) - this%thtr(icell))
      if (thetadif > DEM6) then
        call this%wave_shift(this, icell, icell, -1, &
                             this%nwavst(icell) + 1, 2, -1)
        if (this%uzdpst(2, icell) < DEM30) &
          this%uzdpst(2, icell) = (this%ntrail(icell) + DTWO) * DEM6
        if (this%uzthst(2, icell) > this%thtr(icell)) then
          this%uzspst(2, icell) = this%uzflst(2, icell) / &
                                  (this%uzthst(2, icell) - this%thtr(icell))
        else
          this%uzspst(2, icell) = DZERO
        end if
        this%uzthst(1, icell) = this%thtr(icell)
        this%uzflst(1, icell) = DZERO
        this%uzspst(1, icell) = DZERO
        this%uzdpst(1, icell) = thick
        this%nwavst(icell) = this%nwavst(icell) + 1
        if (this%nwavst(icell) >= this%nwav(icell)) then
          ! -- too many waves error
          ierr = 1
          return
        end if
      else
        this%uzdpst(1, icell) = thick
      end if
    end if
    thetab = this%uzthst(1, icell)
    fluxb = this%uzflst(1, icell)
    this%totflux(icell) = DZERO
    itester = 0
    ffcheck = (this%surflux(icell) - this%uzflst(this%nwavst(icell), icell))
    !
    ! -- increase new waves in infiltration changes
    if (ffcheck > feps2 .OR. ffcheck < -feps2) then
      this%nwavst(icell) = this%nwavst(icell) + 1
      if (this%nwavst(icell) >= this%nwav(icell)) then
        !
        ! -- too many waves error
        ierr = 1
        return
      end if
    else if (this%nwavst(icell) == 1) then
      itester = 1
    end if
    if (this%nwavst(icell) > 1) then
      if (ffcheck < -feps2) then
        call this%trailwav(icell, ierr)
        if (ierr > 0) return
        itrailflg = 1
      end if
      call this%leadwav(time, itester, itrailflg, thetab, fluxb, ffcheck, &
                        feps2, delt, icell)
    end if
    if (itester == 1) then
      this%totflux(icell) = this%totflux(icell) + &
                            (delt - time) * this%uzflst(1, icell)
      time = DZERO
      itester = 0
    end if
    !
    ! -- simulate et
    if (ietflag > 0) call this%uzet(icell, delt, ietflag, ierr)
    if (ierr > 0) return
  end subroutine uzflow

  !> @brief Calculate unit specific tolerances
  !<
  subroutine factors(feps1, feps2)
    ! -- dummy
    real(DP), intent(out) :: feps1
    real(DP), intent(out) :: feps2
    real(DP) :: factor1
    real(DP) :: factor2
    !
    ! calculate constants for uzflow
    factor1 = DONE
    factor2 = DONE
    feps1 = DEM9
    feps2 = DEM9
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
    feps1 = feps1 * factor1 * factor2
    feps2 = feps2 * factor1 * factor2
  end subroutine factors

  !> @brief Create and set trail waves
  !<
  subroutine trailwav(this, icell, ierr)
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    integer(I4B), intent(inout) :: ierr
    ! -- local
    real(DP) :: smoist, smoistinc, ftrail, eps_m1
    real(DP) :: thtsrinv
    real(DP) :: flux1, flux2, theta1, theta2
    real(DP) :: fnuminc
    integer(I4B) :: j, jj, jk, nwavstm1
    !
    ! -- initialize
    eps_m1 = dble(this%eps(icell)) - DONE
    thtsrinv = DONE / (this%thts(icell) - this%thtr(icell))
    nwavstm1 = this%nwavst(icell) - 1
    !
    ! -- initialize trailwaves
    smoist = (((this%surflux(icell) / this%vks(icell))** &
               (DONE / this%eps(icell))) * &
              (this%thts(icell) - this%thtr(icell))) + this%thtr(icell)
    if (this%uzthst(nwavstm1, icell) - smoist > DEM9) then
      fnuminc = DZERO
      do jk = 1, this%ntrail(icell)
        fnuminc = fnuminc + float(jk)
      end do
      smoistinc = (this%uzthst(nwavstm1, icell) - smoist) / (fnuminc - DONE)
      jj = this%ntrail(icell)
      ftrail = dble(this%ntrail(icell)) + DONE
      do j = this%nwavst(icell), this%nwavst(icell) + this%ntrail(icell) - 1
        if (j > this%nwav(icell)) then
          ! -- too many waves error
          ierr = 1
          return
        end if
        if (j > this%nwavst(icell)) then
          this%uzthst(j, icell) = this%uzthst(j - 1, icell) &
                                  - ((ftrail - float(jj)) * smoistinc)
        else
          this%uzthst(j, icell) = this%uzthst(j - 1, icell) - DEM9
        end if
        jj = jj - 1
        if (this%uzthst(j, icell) <= this%thtr(icell) + DEM9) &
          this%uzthst(j, icell) = this%thtr(icell) + DEM9
        this%uzflst(j, icell) = &
          this%vks(icell) * (((this%uzthst(j, icell) - this%thtr(icell)) * &
                              thtsrinv)**this%eps(icell))
        theta2 = this%uzthst(j - 1, icell)
        flux2 = this%uzflst(j - 1, icell)
        flux1 = this%uzflst(j, icell)
        theta1 = this%uzthst(j, icell)
        this%uzspst(j, icell) = leadspeed(theta1, theta2, flux1, flux2, &
                                          this%thts(icell), this%thtr(icell), &
                                          this%eps(icell), this%vks(icell))
        this%uzdpst(j, icell) = DZERO
        if (j == this%nwavst(icell)) then
          this%uzdpst(j, icell) = this%uzdpst(j, icell) + &
                                  (this%ntrail(icell) + 1) * DEM9
        else
          this%uzdpst(j, icell) = this%uzdpst(j - 1, icell) - DEM9
        end if
      end do
      this%nwavst(icell) = this%nwavst(icell) + this%ntrail(icell) - 1
      if (this%nwavst(icell) >= this%nwav(icell)) then
        ! -- too many waves error
        ierr = 1
        return
      end if
    else
      this%uzdpst(this%nwavst, icell) = DZERO
      this%uzflst(this%nwavst, icell) = &
        this%vks(icell) * (((this%uzthst(this%nwavst, icell) - &
                             this%thtr(icell)) * thtsrinv)**this%eps(icell))
      this%uzthst(this%nwavst, icell) = smoist
      theta2 = this%uzthst(this%nwavst(icell) - 1, icell)
      flux2 = this%uzflst(this%nwavst(icell) - 1, icell)
      flux1 = this%uzflst(this%nwavst(icell), icell)
      theta1 = this%uzthst(this%nwavst(icell), icell)
      this%uzspst(this%nwavst(icell), icell) = &
        leadspeed(theta1, theta2, flux1, flux2, this%thts(icell), &
                  this%thtr(icell), this%eps(icell), this%vks(icell))
    end if
  end subroutine trailwav

  !> @brief Create a lead wave and route over time step
  !<
  subroutine leadwav(this, time, itester, itrailflg, thetab, fluxb, &
                     ffcheck, feps2, delt, icell)
    ! -- dummy
    class(UzfCellGroupType) :: this
    real(DP), intent(inout) :: thetab
    real(DP), intent(inout) :: fluxb
    real(DP), intent(in) :: feps2
    real(DP), intent(inout) :: time
    integer(I4B), intent(inout) :: itester
    integer(I4B), intent(inout) :: itrailflg
    real(DP), intent(inout) :: ffcheck
    real(DP), intent(in) :: delt
    integer(I4B), intent(in) :: icell
    ! -- local
    real(DP) :: bottomtime, shortest, fcheck
    real(DP) :: eps_m1, timenew, bottom, timedt
    real(DP) :: thtsrinv, diff, fluxhld2
    real(DP) :: flux1, flux2, theta1, theta2, ftest
    real(DP), allocatable, dimension(:) :: checktime
    integer(I4B) :: iflx, iremove, j, l
    integer(I4B) :: nwavp1, jshort
    integer(I4B), allocatable, dimension(:) :: more
    !
    allocate (checktime(this%nwavst(icell)))
    allocate (more(this%nwavst(icell)))
    ftest = DZERO
    eps_m1 = dble(this%eps(icell)) - DONE
    thtsrinv = DONE / (this%thts(icell) - this%thtr(icell))
    !
    ! -- initialize new wave
    if (itrailflg == 0) then
      if (ffcheck > feps2) then
        this%uzflst(this%nwavst(icell), icell) = this%surflux(icell)
        if (this%uzflst(this%nwavst(icell), icell) < DEM30) &
          this%uzflst(this%nwavst(icell), icell) = DZERO
        this%uzthst(this%nwavst(icell), icell) = &
          (((this%uzflst(this%nwavst(icell), icell) / this%vks(icell))** &
            (DONE / this%eps(icell))) * (this%thts(icell) - this%thtr(icell))) &
          + this%thtr(icell)
        theta2 = this%uzthst(this%nwavst(icell), icell)
        flux2 = this%uzflst(this%nwavst(icell), icell)
        flux1 = this%uzflst(this%nwavst(icell) - 1, icell)
        theta1 = this%uzthst(this%nwavst(icell) - 1, icell)
        this%uzspst(this%nwavst(icell), icell) = &
          leadspeed(theta1, theta2, flux1, flux2, this%thts(icell), &
                    this%thtr(icell), this%eps(icell), this%vks(icell))
        this%uzdpst(this%nwavst(icell), icell) = DZERO
      end if
    end if
    !
    ! -- route all waves and interception of waves over times step
    diff = DONE
    timedt = DZERO
    iflx = 0
    fluxhld2 = this%uzflst(1, icell)
    if (this%nwavst(icell) == 0) itester = 1
    if (itester /= 1) then
      do while (diff > DEM6)
        nwavp1 = this%nwavst(icell) + 1
        timedt = delt - Time
        do j = 1, this%nwavst(icell)
          checktime(j) = DEP20
          more(j) = 0
        end do
        shortest = timedt
        if (this%nwavst(icell) > 2) then
          j = 2
          !
          ! -- calculate time until wave overtakes wave ahead
          nwavp1 = this%nwavst(icell) + 1
          do while (j < nwavp1)
            ftest = this%uzspst(j - 1, icell) - this%uzspst(j, icell)
            if (abs(ftest) > DEM30) then
              checktime(j) = (this%uzdpst(j, icell) - &
                              this%uzdpst(j - 1, icell)) / ftest
              if (checktime(j) < DEM30) checktime(j) = DEP20
            end if
            j = j + 1
          end do
        end if
        !
        ! - calc time until wave reaches bottom of cell
        bottomtime = DEP20
        if (this%nwavst(icell) > 1) then
          if (this%uzspst(2, icell) > DZERO) then
            bottom = this%uzspst(2, icell)
            if (bottom < DEM15) bottom = DEM15
            bottomtime = (this%uzdpst(1, icell) - this%uzdpst(2, icell)) / bottom
            if (bottomtime < DZERO) bottomtime = DEM12
          end if
        end if
        !
        ! -- calc time for wave interception
        jshort = 0
        do j = this%nwavst(icell), 3, -1
          if (shortest - checktime(j) > -DEM9) then
            more(j) = 1
            jshort = j
            shortest = checktime(j)
          end if
        end do
        do j = 3, this%nwavst(icell)
          if (shortest - checktime(j) < DEM9) then
            if (j /= jshort) more(j) = 0
          end if
        end do
        !
        ! -- what happens first, waves hits bottom or interception
        iremove = 0
        timenew = Time
        fcheck = (Time + shortest) - delt
        if (shortest < DEM7) fcheck = -DONE
        if (bottomtime < shortest .AND. Time + bottomtime < delt) then
          j = 2
          do while (j < nwavp1)
            !
            ! -- route waves
            this%uzdpst(j, icell) = this%uzdpst(j, icell) + &
                                    this%uzspst(j, icell) * bottomtime
            j = j + 1
          end do
          fluxb = this%uzflst(2, icell)
          thetab = this%uzthst(2, icell)
          iflx = 1
          call this%wave_shift(this, icell, icell, 1, 1, &
                               this%nwavst(icell) - 1, 1)
          iremove = 1
          timenew = time + bottomtime
          this%uzspst(1, icell) = DZERO
          !
          ! -- do waves intercept before end of time step
        else if (fcheck < DZERO .AND. this%nwavst(icell) > 2) then
          j = 2
          do while (j < nwavp1)
            this%uzdpst(j, icell) = this%uzdpst(j, icell) + &
                                    this%uzspst(j, icell) * shortest
            j = j + 1
          end do
          !
          ! -- combine waves that intercept, remove a wave
          j = 3
          l = j
          do while (j < this%nwavst(icell) + 1)
            if (more(j) == 1) then
              l = j
              theta2 = this%uzthst(j, icell)
              flux2 = this%uzflst(j, icell)
              if (j == 3) then
                flux1 = fluxb
                theta1 = thetab
              else
                flux1 = this%uzflst(j - 2, icell)
                theta1 = this%uzthst(j - 2, icell)
              end if
              this%uzspst(j, icell) = leadspeed(theta1, theta2, flux1, flux2, &
                                                this%thts(icell), &
                                                this%thtr(icell), &
                                                this%eps(icell), this%vks(icell))
              !
              ! -- update waves
              call this%wave_shift(this, icell, icell, 1, l - 1, &
                                   this%nwavst(icell) - 1, 1)
              l = this%nwavst(icell) + 1
              iremove = iremove + 1
            end if
            j = j + 1
          end do
          timenew = timenew + shortest
          !
          ! -- calc. total flux to bottom during remaining time in step
        else
          j = 2
          do while (j < nwavp1)
            this%uzdpst(j, icell) = this%uzdpst(j, icell) + &
                                    this%uzspst(j, icell) * timedt
            j = j + 1
          end do
          timenew = delt
        end if
        this%totflux(icell) = this%totflux(icell) + fluxhld2 * (timenew - time)
        if (iflx == 1) then
          fluxhld2 = this%uzflst(1, icell)
          iflx = 0
        end if
        !
        ! -- remove dead waves
        this%nwavst(icell) = this%nwavst(icell) - iremove
        time = timenew
        diff = delt - Time
        if (this%nwavst(icell) == 1) then
          itester = 1
          exit
        end if
      end do
    end if
    deallocate (checktime)
    deallocate (more)
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
    j = this%nwavst(icell) + 1
    k = this%nwavst(icell)
    nwavm1 = k - 1
    if (d1 > this%uzdpst(1, icell)) d1 = this%uzdpst(1, icell)
    !
    ! -- find deepest wave above depth d1, counter held as j
    do while (k > 0)
      if (this%uzdpst(k, icell) - d1 < -DEM30) j = k
      k = k - 1
    end do
    if (j > this%nwavst(icell)) then
      fm = fm + (this%uzthst(this%nwavst(icell), icell) - this%thtr(icell)) * d1
    elseif (this%nwavst(icell) > 1) then
      if (j > 1) then
        fm = fm + (this%uzthst(j - 1, icell) - this%thtr(icell)) &
             * (d1 - this%uzdpst(j, icell))
      end if
      do jj = j, nwavm1
        fm = fm + (this%uzthst(jj, icell) - this%thtr(icell)) &
             * (this%uzdpst(jj, icell) &
                - this%uzdpst(jj + 1, icell))
      end do
      fm = fm + (this%uzthst(this%nwavst(icell), icell) - this%thtr(icell)) &
           * (this%uzdpst(this%nwavst(icell), icell))
    else
      fm = fm + (this%uzthst(1, icell) - this%thtr(icell)) * d1
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
    real(DP) :: thick, thtsrinv
    integer(I4B) :: nwavhld, k, j
    !
    bot = this%watab(icell)
    top = this%celtop(icell)
    thick = top - bot
    nwavhld = this%nwavst(icell)
    if (itest == 1) then
      this%uzflst(1, icell) = DZERO
      this%uzthst(1, icell) = this%thtr(icell)
      return
    end if
    if (iss == 1) then
      if (this%thts(icell) - this%thtr(icell) < DEM7) then
        thtsrinv = DONE / DEM7
      else
        thtsrinv = DONE / (this%thts(icell) - this%thtr(icell))
      end if
      this%totflux(icell) = this%surflux(icell) * delt
      this%watabold(icell) = this%watab(icell)
      this%uzthst(1, icell) = this%thti(icell)
      this%uzflst(1, icell) = &
        this%vks(icell) * (((this%uzthst(1, icell) - this%thtr(icell)) &
                            * thtsrinv)**this%eps(icell))
      this%uzdpst(1, icell) = thick
      this%uzspst(1, icell) = thick
      this%nwavst(icell) = 1
    else
      !
      ! -- water table rises through waves
      if (this%watab(icell) - this%watabold(icell) > DEM30) then
        depthsave = this%uzdpst(1, icell)
        j = 0
        k = this%nwavst(icell)
        do while (k > 0)
          if (this%uzdpst(k, icell) - thick < -DEM30) j = k
          k = k - 1
        end do
        this%uzdpst(1, icell) = thick
        if (j > 1) then
          this%uzspst(1, icell) = DZERO
          this%nwavst(icell) = this%nwavst(icell) - j + 2
          this%uzthst(1, icell) = this%uzthst(j - 1, icell)
          this%uzflst(1, icell) = this%uzflst(j - 1, icell)
          if (j > 2) call this%wave_shift(this, icell, icell, j - 2, 2, &
                                          nwavhld - (j - 2), 1)
        elseif (j == 0) then
          this%uzspst(1, icell) = DZERO
          this%uzthst(1, icell) = this%uzthst(this%nwavst(icell), icell)
          this%uzflst(1, icell) = this%uzflst(this%nwavst(icell), icell)
          this%nwavst(icell) = 1
        end if
      end if
      !
      ! -- calculate new unsat. storage
      if (thick <= DZERO) then
        this%uzspst(1, icell) = DZERO
        this%nwavst(icell) = 1
        this%uzthst(1, icell) = this%thtr(icell)
        this%uzflst(1, icell) = DZERO
      end if
      this%watabold(icell) = this%watab(icell)
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
    type(UzfCellGroupType) :: uzfktemp
    real(DP) :: diff
    real(DP) :: thetaout
    real(DP) :: fm
    real(DP) :: st
    real(DP) :: thtsrinv
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
    this%etact(icell) = DZERO
    if (this%extdpuz(icell) < DEM7) return
    petsub = this%rootact(icell) * this%pet(icell) * &
             this%extdpuz(icell) / this%extdp(icell)
    thetaout = delt * petsub / this%extdp(icell)
    if (ietflag == 1) thetaout = delt * this%pet(icell) / this%extdp(icell)
    if (thetaout < DEM10) return
    depth = this%uzdpst(1, icell)
    st = this%unsat_stor(icell, depth)
    if (st < DEM4) return
    !
    ! -- allocate temporary wave storage.
    nwv = this%nwavst(icell)
    itest = 0
    call uzfktemp%init(1, nwv)
    !
    ! store original wave characteristics
    call uzfktemp%wave_shift(this, 1, icell, 0, 1, nwv, 1)
    factor = DONE
    this%etact(icell) = DZERO
    if (this%thts(icell) - this%thtr(icell) < DEM7) then
      thtsrinv = 1.0 / DEM7
    else
      thtsrinv = DONE / (this%thts(icell) - this%thtr(icell))
    end if
    epsfksthts = this%eps(icell) * this%vks(icell) * thtsrinv
    this%etact(icell) = DZERO
    fmp = DZERO
    extwc1 = this%extwc(icell) - this%thtr(icell)
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
      if (this%nwavst(icell) == 1 .AND. &
          this%uzdpst(1, icell) <= this%extdpuz(icell)) then
        if (ietflag == 2) then
          tho = this%uzthst(1, icell)
          fktho = this%uzflst(1, icell)
          hcap = this%caph(icell, tho)
          thetaout = this%rate_et_z(icell, factor, fktho, hcap)
        end if
        if ((this%uzthst(1, icell) - thetaout) > this%thtr(icell) + extwc1) then
          this%uzthst(1, icell) = this%uzthst(1, icell) - thetaout
          this%uzflst(1, icell) = &
            this%vks(icell) * (((this%uzthst(1, icell) - &
                                 this%thtr(icell)) * thtsrinv)**this%eps(icell))
        else if (this%uzthst(1, icell) > this%thtr(icell) + extwc1) then
          this%uzthst(1, icell) = this%thtr(icell) + extwc1
          this%uzflst(1, icell) = &
            this%vks(icell) * (((this%uzthst(1, icell) - &
                                 this%thtr(icell)) * thtsrinv)**this%eps(icell))
        end if
        !
        ! -- all waves shallower than extinction depth
      else if (this%nwavst(icell) > 1 .AND. &
               this%uzdpst(this%nwavst(icell), icell) > this%extdpuz(icell)) then
        if (ietflag == 2) then
          tho = this%uzthst(this%nwavst(icell), icell)
          fktho = this%uzflst(this%nwavst(icell), icell)
          hcap = this%caph(icell, tho)
          thetaout = this%rate_et_z(icell, factor, fktho, hcap)
        end if
        if (this%uzthst(this%nwavst(icell), icell) - thetaout > &
            this%thtr(icell) + extwc1) then
          this%uzthst(this%nwavst(icell) + 1, icell) = &
            this%uzthst(this%nwavst(icell), icell) - thetaout
          numadd = 1
        else if (this%uzthst(this%nwavst(icell), icell) > &
                 this%thtr(icell) + extwc1) then
          this%uzthst(this%nwavst(icell) + 1, icell) = this%thtr(icell) + extwc1
          numadd = 1
        end if
        if (numadd == 1) then
          this%uzflst(this%nwavst(icell) + 1, icell) = &
            this%vks(icell) * &
            (((this%uzthst(this%nwavst(icell) + 1, icell) - &
               this%thtr(icell)) * thtsrinv)**this%eps(icell))
          theta2 = this%uzthst(this%nwavst(icell) + 1, icell)
          flux2 = this%uzflst(this%nwavst(icell) + 1, icell)
          flux1 = this%uzflst(this%nwavst(icell), icell)
          theta1 = this%uzthst(this%nwavst(icell), icell)
          this%uzspst(this%nwavst(icell) + 1, icell) = &
            leadspeed(theta1, theta2, flux1, flux2, this%thts(icell), &
                      this%thtr(icell), this%eps(icell), this%vks(icell))
          this%uzdpst(this%nwavst(icell) + 1, icell) = this%extdpuz(icell)
          this%nwavst(icell) = this%nwavst(icell) + 1
          if (this%nwavst(icell) > this%nwav(icell)) then
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
      else if (this%nwavst(icell) == 1) then
        if (ietflag == 2) then
          tho = this%uzthst(1, icell)
          fktho = this%uzflst(1, icell)
          hcap = this%caph(icell, tho)
          thetaout = this%rate_et_z(icell, factor, fktho, hcap)
        end if
        if ((this%uzthst(1, icell) - thetaout) > this%thtr(icell) + extwc1) then
          if (thetaout > DEM30) then
            this%uzthst(2, icell) = this%uzthst(1, icell) - thetaout
            this%uzflst(2, icell) = &
              this%vks(icell) * (((this%uzthst(2, icell) - this%thtr(icell)) * &
                                  thtsrinv)**this%eps(icell))
            this%uzdpst(2, icell) = this%extdpuz(icell)
            theta2 = this%uzthst(2, icell)
            flux2 = this%uzflst(2, icell)
            flux1 = this%uzflst(1, icell)
            theta1 = this%uzthst(1, icell)
            this%uzspst(2, icell) = &
              leadspeed(theta1, theta2, flux1, flux2, this%thts(icell), &
                        this%thtr(icell), this%eps(icell), this%vks(icell))
            this%nwavst(icell) = this%nwavst(icell) + 1
            if (this%nwavst(icell) > this%nwav(icell)) then
              !
              ! -- too many waves error
              ierr = 1
              goto 500
            end if
          end if
        else if (this%uzthst(1, icell) > this%thtr(icell) + extwc1) then
          if (thetaout > DEM30) then
            this%uzthst(2, icell) = this%thtr(icell) + extwc1
            this%uzflst(2, icell) = &
              this%vks(icell) * (((this%uzthst(2, icell) - &
                                   this%thtr(icell)) * thtsrinv)**this%eps(icell))
            this%uzdpst(2, icell) = this%extdpuz(icell)
            theta2 = this%uzthst(2, icell)
            flux2 = this%uzflst(2, icell)
            flux1 = this%uzflst(1, icell)
            theta1 = this%uzthst(1, icell)
            this%uzspst(2, icell) = &
              leadspeed(theta1, theta2, flux1, flux2, this%thts(icell), &
                        this%thtr(icell), this%eps(icell), this%vks(icell))
            this%nwavst(icell) = this%nwavst(icell) + 1
            if (this%nwavst(icell) > this%nwav(icell)) then
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
        if (this%uzdpst(1, icell) - this%extdpuz(icell) > DEM7) then
          j = 2
          jk = 0
          !
          ! -- locate extinction depth between waves
          do while (jk == 0)
            diff = this%uzdpst(j, icell) - this%extdpuz(icell)
            if (diff > dzero) then
              j = j + 1
            else
              jk = 1
            end if
          end do
          kk = j
          if (this%uzthst(j, icell) > this%thtr(icell) + extwc1) then
            !
            ! -- create a wave at extinction depth
            if (abs(diff) > DEM5) then
              call this%wave_shift(this, icell, icell, -1, &
                                   this%nwavst(icell) + 1, j, -1)
              this%uzdpst(j, icell) = this%extdpuz(icell)
              this%nwavst(icell) = this%nwavst(icell) + 1
              if (this%nwavst(icell) > this%nwav(icell)) then
                !
                ! -- too many waves error
                ierr = 1
                goto 500
              end if
            end if
            kk = j
          else
            jhold = this%nwavst(icell)
            i = j + 1
            do while (i < this%nwavst(icell))
              if (this%uzthst(i, icell) > this%thtr(icell) + extwc1) then
                jhold = i
                i = this%nwavst(icell) + 1
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
        do while (kk <= this%nwavst(icell))
          if (ietflag == 2) then
            tho = this%uzthst(kk, icell)
            fktho = this%uzflst(kk, icell)
            hcap = this%caph(icell, tho)
            thetaout = this%rate_et_z(icell, factor, fktho, hcap)
          end if
          if (this%uzthst(kk, icell) > this%thtr(icell) + extwc1) then
            if (this%uzthst(kk, icell) - thetaout > &
                this%thtr(icell) + extwc1) then
              this%uzthst(kk, icell) = this%uzthst(kk, icell) - thetaout
            else if (this%uzthst(kk, icell) > this%thtr(icell) + extwc1) then
              this%uzthst(kk, icell) = this%thtr(icell) + extwc1
            end if
            if (kk == 1) then
              this%uzflst(kk, icell) = &
                this%vks(icell) * &
                (((this%uzthst(kk, icell) - &
                   this%thtr(icell)) * thtsrinv)**this%eps(icell))
            end if
            if (kk > 1) then
              flux1 = &
                this%vks(icell) * ((this%uzthst(kk - 1, icell) - &
                                    this%thtr(icell)) * thtsrinv)**this%eps(icell)
              flux2 = &
                this%vks(icell) * ((this%uzthst(kk, icell) - &
                                    this%thtr(icell)) * thtsrinv)**this%eps(icell)
              this%uzflst(kk, icell) = flux2
              theta2 = this%uzthst(kk, icell)
              theta1 = this%uzthst(kk - 1, icell)
              this%uzspst(kk, icell) = leadspeed(theta1, theta2, flux1, flux2, &
                                                 this%thts(icell), &
                                                 this%thtr(icell), &
                                                 this%eps(icell), this%vks(icell))
            end if
          end if
          kk = kk + 1
        end do
      end if
      !
      ! -- calculate aet
      kj = 1
      do while (kj <= this%nwavst(icell) - 1)
        if (abs(this%uzthst(kj, icell) - this%uzthst(kj + 1, icell)) < DEM6) then
          call this%wave_shift(this, icell, icell, 1, kj + 1, &
                               this%nwavst(icell) - 1, 1)
          kj = kj - 1
          this%nwavst(icell) = this%nwavst(icell) - 1
        end if
        kj = kj + 1
      end do
      depth = this%uzdpst(1, icell)
      fm = this%unsat_stor(icell, depth)
      this%etact(icell) = st - fm
      fm = this%etact(icell) / delt
      if (this%etact(icell) < dzero) then
        call this%wave_shift(uzfktemp, icell, 1, 0, 1, nwv, 1)
        this%nwavst(icell) = nwv
        this%etact(icell) = DZERO
      elseif (petsub - fm < -DEM15 .AND. ietflag == 2) then
        !
        ! -- aet greater than pet, reset and try again
        call this%wave_shift(uzfktemp, icell, 1, 0, 1, nwv, 1)
        this%nwavst(icell) = nwv
        this%etact(icell) = DZERO
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
    !
    ! -- deallocate temporary worker
    call uzfktemp%dealloc()
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
    star = (tho - this%thtr(icell)) / (this%thts(icell) - this%thtr(icell))
    if (star < DEM15) star = DEM15
    lambda = DTWO / (this%eps(icell) - DTHREE)
    if (star > DEM15) then
      if (tho - this%thts(icell) < DEM15) then
        caph = this%ha(icell) * star**(-DONE / lambda)
      else
        caph = DZERO
      end if
    end if
  end function caph

  !> @brief Calculate capillary pressure-based uz et
  function rate_et_z(this, icell, factor, fktho, h)
    ! -- Return
    real(DP) :: rate_et_z
    ! -- dummy
    class(UzfCellGroupType) :: this
    integer(I4B), intent(in) :: icell
    real(DP), intent(in) :: factor, fktho, h
    !
    rate_et_z = factor * fktho * (h - this%hroot(icell))
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
    if (this%watab(icell) < this%celtop(icell)) then
      if (this%celtop(icell) - depth > this%watab(icell)) then
        d1 = depth - DEM3
        d2 = depth + DEM3
        f1 = this%unsat_stor(icell, d1)
        f2 = this%unsat_stor(icell, d2)
        theta_at_depth = this%thtr(icell) + (f2 - f1) / (d2 - d1)
      else
        theta_at_depth = this%thts(icell)
      end if
    else
      theta_at_depth = this%thts(icell)
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
    hgwf = this%watab(icell)
    top = this%celtop(icell)
    bot = this%celbot(icell)
    thk = top - max(bot, hgwf)
    if (thk > DZERO) then
      theta_r = this%thtr(icell)
      d = thk
      fm = this%unsat_stor(icell, d)
      watercontent = fm / thk
      watercontent = watercontent + theta_r
    else
      watercontent = DZERO
    end if
  end function get_wcnew

end module UzfCellGroupModule
