! -- Lake Energy Transport Module
! -- todo: what to do about reactions in lake?  Decay?
! -- todo: save the lke temperatures into the lak aux variable?
!
! LAK flows (lakbudptr)     index var     LKT term              Transport Type
!---------------------------------------------------------------------------------

! -- terms from LAK that will be handled by parent APT Package
! FLOW-JA-FACE              idxbudfjf     FLOW-JA-FACE          cv2cv
! GWF (aux FLOW-AREA)       idxbudgwf     GWF                   cv2gwf
! STORAGE (aux VOLUME)      idxbudsto     none                  used for cv volumes
! FROM-MVR                  idxbudfmvr    FROM-MVR              q * cext = this%qfrommvr(:)
! TO-MVR                    idxbudtmvr    TO-MVR                q * cfeat

! -- LAK terms
! RAINFALL                  idxbudrain    RAINFALL              q * crain
! EVAPORATION               idxbudevap    EVAPORATION           cfeat<cevap: q*cfeat, else: q*cevap
! RUNOFF                    idxbudroff    RUNOFF                q * croff
! EXT-INFLOW                idxbudiflw    EXT-INFLOW            q * ciflw
! WITHDRAWAL                idxbudwdrl    WITHDRAWAL            q * cfeat
! EXT-OUTFLOW               idxbudoutf    EXT-OUTFLOW           q * cfeat
! LAKEBED-COND              idxbudlbcd    LAKEBED-COND          ! kluge note: expression for this

! -- terms from a flow file that should be skipped
! CONSTANT                  none          none                  none
! AUXILIARY                 none          none                  none

! -- terms that are written to the transport budget file
! none                      none          STORAGE (aux MASS)    dM/dt
! none                      none          AUXILIARY             none
! none                      none          CONSTANT              accumulate
!
!
module GweLkeModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, LINELENGTH, LENBOUNDNAME, DEP20
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors
  use BndModule, only: BndType, GetBndFromList
  use TspFmiModule, only: TspFmiType
  use LakModule, only: LakType
  use ObserveModule, only: ObserveType
  use TspAptModule, only: TspAptType, apt_process_obsID, &
                          apt_process_obsID12
  use GweInputDataModule, only: GweInputDataType
  use MatrixBaseModule

  implicit none

  public lke_create

  character(len=*), parameter :: ftype = 'LKE'
  character(len=*), parameter :: flowtype = 'LAK'
  character(len=16) :: text = '             LKE'

  type, extends(TspAptType) :: GweLkeType

    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in mst

    integer(I4B), pointer :: idxbudrain => null() ! index of rainfall terms in flowbudptr
    integer(I4B), pointer :: idxbudevap => null() ! index of evaporation terms in flowbudptr
    integer(I4B), pointer :: idxbudroff => null() ! index of runoff terms in flowbudptr
    integer(I4B), pointer :: idxbudiflw => null() ! index of inflow terms in flowbudptr
    integer(I4B), pointer :: idxbudwdrl => null() ! index of withdrawal terms in flowbudptr
    integer(I4B), pointer :: idxbudoutf => null() ! index of outflow terms in flowbudptr
    integer(I4B), pointer :: idxbudlbcd => null() ! index of lakebed conduction terms in flowbudptr

    real(DP), dimension(:), pointer, contiguous :: temprain => null() ! rainfall temperature
    real(DP), dimension(:), pointer, contiguous :: tempevap => null() ! evaporation temperature
    real(DP), dimension(:), pointer, contiguous :: temproff => null() ! runoff temperature
    real(DP), dimension(:), pointer, contiguous :: tempiflw => null() ! inflow temperature
    real(DP), dimension(:), pointer, contiguous :: ktf => null() !< thermal conductivity between the lke and groundwater cell
    real(DP), dimension(:), pointer, contiguous :: rfeatthk => null() !< thickness of lakebed material through which thermal conduction occurs

  contains

    procedure :: bnd_da => lke_da
    procedure :: allocate_scalars
    procedure :: apt_allocate_arrays => lke_allocate_arrays
    procedure :: find_apt_package => find_lke_package
    procedure :: pak_fc_expanded => lke_fc_expanded
    procedure :: pak_solve => lke_solve
    procedure :: pak_get_nbudterms => lke_get_nbudterms
    procedure :: pak_setup_budobj => lke_setup_budobj
    procedure :: pak_fill_budobj => lke_fill_budobj
    procedure :: lke_rain_term
    procedure :: lke_evap_term
    procedure :: lke_roff_term
    procedure :: lke_iflw_term
    procedure :: lke_wdrl_term
    procedure :: lke_outf_term
    procedure :: pak_df_obs => lke_df_obs
    procedure :: pak_rp_obs => lke_rp_obs
    procedure :: pak_bd_obs => lke_bd_obs
    procedure :: pak_set_stressperiod => lke_set_stressperiod
    procedure :: apt_read_cvs => lke_read_cvs

  end type GweLkeType

contains

  !> @brief Create a new lke package
  !<
  subroutine lke_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        fmi, eqnsclfac, gwecommon, dvt, dvu, dvua)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    type(TspFmiType), pointer :: fmi
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    type(GweInputDataType), intent(in), target :: gwecommon !< shared data container for use by multiple GWE packages
    character(len=*), intent(in) :: dvt !< For GWE, set to "TEMPERATURE" in TspAptType
    character(len=*), intent(in) :: dvu !< For GWE, set to "energy" in TspAptType
    character(len=*), intent(in) :: dvua !< For GWE, set to "E" in TspAptType
    ! -- local
    type(GweLkeType), pointer :: lkeobj
    !
    ! -- Allocate the object and assign values to object variables
    allocate (lkeobj)
    packobj => lkeobj
    !
    ! -- Create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- Allocate scalars
    call lkeobj%allocate_scalars()
    !
    ! -- Initialize package
    call packobj%pack_initialize()
    !
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1
    !
    ! -- Store pointer to flow model interface.  When the GwfGwe exchange is
    !    created, it sets fmi%bndlist so that the GWE model has access to all
    !    the flow packages
    lkeobj%fmi => fmi
    !
    ! -- Store pointer to governing equation scale factor
    lkeobj%eqnsclfac => eqnsclfac
    !
    ! -- Store pointer to shared data module for accessing cpw, rhow
    !    for the budget calculations, and for accessing the latent heat of
    !    vaporization for evaporative cooling.
    lkeobj%gwecommon => gwecommon
    !
    ! -- Set labels that will be used in generalized APT class
    lkeobj%depvartype = dvt
    lkeobj%depvarunit = dvu
    lkeobj%depvarunitabbrev = dvua
  end subroutine lke_create

  !> @brief Find corresponding lke package
  !<
  subroutine find_lke_package(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweLkeType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg
    class(BndType), pointer :: packobj
    integer(I4B) :: ip, icount
    integer(I4B) :: nbudterm
    logical :: found
    !
    ! -- Initialize found to false, and error later if flow package cannot
    !    be found
    found = .false.
    !
    ! -- If user is specifying flows in a binary budget file, then set up
    !    the budget file reader, otherwise set a pointer to the flow package
    !    budobj
    if (this%fmi%flows_from_file) then
      call this%fmi%set_aptbudobj_pointer(this%flowpackagename, this%flowbudptr)
      if (associated(this%flowbudptr)) found = .true.
      !
    else
      if (associated(this%fmi%gwfbndlist)) then
        ! -- Look through gwfbndlist for a flow package with the same name as
        !    this transport package name
        do ip = 1, this%fmi%gwfbndlist%Count()
          packobj => GetBndFromList(this%fmi%gwfbndlist, ip)
          if (packobj%packName == this%flowpackagename) then
            found = .true.
            !
            ! -- Store BndType pointer to packobj, and then
            !    use the select type to point to the budobj in flow package
            this%flowpackagebnd => packobj
            select type (packobj)
            type is (LakType)
              this%flowbudptr => packobj%budobj
            end select
          end if
          if (found) exit
        end do
      end if
    end if
    !
    ! -- Error if flow package not found
    if (.not. found) then
      write (errmsg, '(a)') 'Could not find flow package with name '&
                            &//trim(adjustl(this%flowpackagename))//'.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Allocate space for idxbudssm, which indicates whether this is a
    !    special budget term or one that is a general source and sink
    nbudterm = this%flowbudptr%nbudterm
    call mem_allocate(this%idxbudssm, nbudterm, 'IDXBUDSSM', this%memoryPath)
    !
    ! -- Process budget terms and identify special budget terms
    write (this%iout, '(/, a, a)') &
      'PROCESSING '//ftype//' INFORMATION FOR ', this%packName
    write (this%iout, '(a)') '  IDENTIFYING FLOW TERMS IN '//flowtype//' PACKAGE'
    write (this%iout, '(a, i0)') &
      '  NUMBER OF '//flowtype//' = ', this%flowbudptr%ncv
    icount = 1
    do ip = 1, this%flowbudptr%nbudterm
      select case (trim(adjustl(this%flowbudptr%budterm(ip)%flowtype)))
      case ('FLOW-JA-FACE')
        this%idxbudfjf = ip
        this%idxbudssm(ip) = 0
      case ('GWF')
        this%idxbudgwf = ip
        this%idxbudlbcd = ip
        this%idxbudssm(ip) = 0
      case ('STORAGE')
        this%idxbudsto = ip
        this%idxbudssm(ip) = 0
      case ('RAINFALL')
        this%idxbudrain = ip
        this%idxbudssm(ip) = 0
      case ('EVAPORATION')
        this%idxbudevap = ip
        this%idxbudssm(ip) = 0
      case ('RUNOFF')
        this%idxbudroff = ip
        this%idxbudssm(ip) = 0
      case ('EXT-INFLOW')
        this%idxbudiflw = ip
        this%idxbudssm(ip) = 0
      case ('WITHDRAWAL')
        this%idxbudwdrl = ip
        this%idxbudssm(ip) = 0
      case ('EXT-OUTFLOW')
        this%idxbudoutf = ip
        this%idxbudssm(ip) = 0
      case ('TO-MVR')
        this%idxbudtmvr = ip
        this%idxbudssm(ip) = 0
      case ('FROM-MVR')
        this%idxbudfmvr = ip
        this%idxbudssm(ip) = 0
      case ('AUXILIARY')
        this%idxbudaux = ip
        this%idxbudssm(ip) = 0
      case default
        !
        ! -- Set idxbudssm equal to a column index for where the temperatures
        !    are stored in the tempbud(nbudssm, ncv) array
        this%idxbudssm(ip) = icount
        icount = icount + 1
      end select
      write (this%iout, '(a, i0, " = ", a,/, a, i0)') &
        '  TERM ', ip, trim(adjustl(this%flowbudptr%budterm(ip)%flowtype)), &
        '   MAX NO. OF ENTRIES = ', this%flowbudptr%budterm(ip)%maxlist
    end do
    write (this%iout, '(a, //)') 'DONE PROCESSING '//ftype//' INFORMATION'
  end subroutine find_lke_package

  !> @brief Add matrix terms related to LKE
  !!
  !! This will be called from TspAptType%apt_fc_expanded()
  !! in order to add matrix terms specifically for LKE
  !<
  subroutine lke_fc_expanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GweLkeType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, n, n1, n2
    integer(I4B) :: iloc
    integer(I4B) :: iposd, iposoffd
    integer(I4B) :: ipossymd, ipossymoffd
    integer(I4B) :: auxpos
    real(DP) :: rrate
    real(DP) :: rhsval
    real(DP) :: hcofval
    real(DP) :: ctherm !< thermal conductance
    real(DP) :: wa !< wetted area
    real(DP) :: ktf !< thermal conductivity of streambed material
    real(DP) :: s !< thickness of conductive streambed material
    !
    ! -- Add rainfall contribution
    if (this%idxbudrain /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrain)%nlist
        call this%lke_rain_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add evaporation contribution
    if (this%idxbudevap /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudevap)%nlist
        call this%lke_evap_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add runoff contribution
    if (this%idxbudroff /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudroff)%nlist
        call this%lke_roff_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add inflow contribution
    if (this%idxbudiflw /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudiflw)%nlist
        call this%lke_iflw_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add withdrawal contribution
    if (this%idxbudwdrl /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudwdrl)%nlist
        call this%lke_wdrl_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add outflow contribution
    if (this%idxbudoutf /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudoutf)%nlist
        call this%lke_outf_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add lakebed conduction contribution
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      !
      ! -- Set n to feature number and process if active feature
      n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      if (this%iboundpak(n) /= 0) then
        !
        ! -- Set acoef and rhs to negative so they are relative to sfe and not gwe
        auxpos = this%flowbudptr%budterm(this%idxbudgwf)%naux
        wa = this%flowbudptr%budterm(this%idxbudgwf)%auxvar(auxpos, j)
        ktf = this%ktf(n)
        s = this%rfeatthk(n)
        ctherm = ktf * wa / s
        !
        ! -- Add to sfe row
        iposd = this%idxdglo(j)
        iposoffd = this%idxoffdglo(j)
        call matrix_sln%add_value_pos(iposd, -ctherm) ! kluge note: make sure the signs on ctherm are correct here and below
        call matrix_sln%add_value_pos(iposoffd, ctherm)
        !
        ! -- Add to gwe row for sfe connection
        ipossymd = this%idxsymdglo(j)
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymd, -ctherm)
        call matrix_sln%add_value_pos(ipossymoffd, ctherm)
      end if
    end do
  end subroutine lke_fc_expanded

  !> @brief Add terms specific to lakes to the explicit lake solve
  !<
  subroutine lke_solve(this)
    ! -- dummy
    class(GweLkeType) :: this
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n1, n2
    real(DP) :: rrate
    !
    ! -- Add rainfall contribution
    if (this%idxbudrain /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrain)%nlist
        call this%lke_rain_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add evaporation contribution
    if (this%idxbudevap /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudevap)%nlist
        call this%lke_evap_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add runoff contribution
    if (this%idxbudroff /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudroff)%nlist
        call this%lke_roff_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add inflow contribution
    if (this%idxbudiflw /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudiflw)%nlist
        call this%lke_iflw_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add withdrawal contribution
    if (this%idxbudwdrl /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudwdrl)%nlist
        call this%lke_wdrl_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add outflow contribution
    if (this%idxbudoutf /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudoutf)%nlist
        call this%lke_outf_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
  end subroutine lke_solve

  !> @brief Function to return the number of budget terms just for this package.
  !!
  !! This overrides a function in the parent class.
  !<
  function lke_get_nbudterms(this) result(nbudterms)
    ! -- dummy
    class(GweLkeType) :: this
    ! -- Return
    integer(I4B) :: nbudterms
    !
    ! -- Number of budget terms is 7
    !    1) rainfall
    !    2) evap
    !    3) runoff
    !    4) ext-inflow
    !    5) withdrawal
    !    6) ext-outflow
    !    7) lakebed-cond
    !
    nbudterms = 7
  end function lke_get_nbudterms

  !> @brief Set up the budget object that stores all the lake flows
  !<
  subroutine lke_setup_budobj(this, idx)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(inout) :: idx
    ! -- local
    integer(I4B) :: n, n1, n2
    integer(I4B) :: maxlist, naux
    real(DP) :: q
    character(len=LENBUDTXT) :: text
    !
    ! -- Addition of heat associated with rainfall directly on the lake surface
    text = '        RAINFALL'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudrain)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Evaporative cooling from lake surface
    text = '     EVAPORATION'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudevap)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Addition of heat associated with runoff that flows to the lake
    text = '          RUNOFF'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudroff)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Addition of heat associated with user-specified inflow to the lake
    text = '      EXT-INFLOW'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudiflw)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Removal of heat associated with user-specified withdrawal from lake
    text = '      WITHDRAWAL'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudwdrl)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Removal of heat associated with outflow from lake that leaves
    !    model domain
    text = '     EXT-OUTFLOW'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudoutf)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Conductive exchange of heat through the wetted lakebed
    text = '    LAKEBED-COND'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudlbcd)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    call this%budobj%budterm(idx)%reset(maxlist)
    q = DZERO
    do n = 1, maxlist
      n1 = this%flowbudptr%budterm(this%idxbudgwf)%id1(n)
      n2 = this%flowbudptr%budterm(this%idxbudgwf)%id2(n)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
    end do
  end subroutine lke_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  subroutine lke_fill_budobj(this, idx, x, flowja, ccratin, ccratout)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(inout) :: idx
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    real(DP), intent(inout) :: ccratin
    real(DP), intent(inout) :: ccratout
    ! -- local
    integer(I4B) :: j, n1, n2
    integer(I4B) :: nlist
    integer(I4B) :: igwfnode
    integer(I4B) :: idiag
    integer(I4B) :: auxpos
    real(DP) :: q
    real(DP) :: ctherm !< thermal conductance
    real(DP) :: wa !< wetted area
    real(DP) :: ktf !< thermal conductivity of streambed material
    real(DP) :: s !< thickness of conductive streambed materia
    !
    ! -- Rain
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudrain)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%lke_rain_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- Evaporation
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudevap)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%lke_evap_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- Runoff
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudroff)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%lke_roff_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- Est-Inflow
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudiflw)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%lke_iflw_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- Withdrawal
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudwdrl)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%lke_wdrl_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- Ext-Outflow
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudoutf)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%lke_outf_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- Lakebed-Cond
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do j = 1, this%flowbudptr%budterm(this%idxbudlbcd)%nlist
      q = DZERO
      n1 = this%flowbudptr%budterm(this%idxbudlbcd)%id1(j)
      if (this%iboundpak(n1) /= 0) then
        igwfnode = this%flowbudptr%budterm(this%idxbudlbcd)%id2(j)
        auxpos = this%flowbudptr%budterm(this%idxbudgwf)%naux ! for now there is only 1 aux variable under 'GWF'
        wa = this%flowbudptr%budterm(this%idxbudgwf)%auxvar(auxpos, j)
        ktf = this%ktf(n1)
        s = this%rfeatthk(n1)
        ctherm = ktf * wa / s
        q = ctherm * (x(igwfnode) - this%xnewpak(n1))
      end if
      call this%budobj%budterm(idx)%update_term(n1, igwfnode, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      if (this%iboundpak(n1) /= 0) then
        ! -- Contribution to gwe cell budget
        this%simvals(n1) = this%simvals(n1) - q
        idiag = this%dis%con%ia(igwfnode)
        flowja(idiag) = flowja(idiag) - q
      end if
    end do
  end subroutine lke_fill_budobj

  !> @brief Allocate scalars specific to the lake energy transport (LKE)
  !! package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweLkeType) :: this
    !
    ! -- Allocate scalars in TspAptType
    call this%TspAptType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%idxbudrain, 'IDXBUDRAIN', this%memoryPath)
    call mem_allocate(this%idxbudevap, 'IDXBUDEVAP', this%memoryPath)
    call mem_allocate(this%idxbudroff, 'IDXBUDROFF', this%memoryPath)
    call mem_allocate(this%idxbudiflw, 'IDXBUDIFLW', this%memoryPath)
    call mem_allocate(this%idxbudwdrl, 'IDXBUDWDRL', this%memoryPath)
    call mem_allocate(this%idxbudoutf, 'IDXBUDOUTF', this%memoryPath)
    call mem_allocate(this%idxbudlbcd, 'IDXBUDLBCD', this%memoryPath)
    !
    ! -- Initialize
    this%idxbudrain = 0
    this%idxbudevap = 0
    this%idxbudroff = 0
    this%idxbudiflw = 0
    this%idxbudwdrl = 0
    this%idxbudoutf = 0
    this%idxbudlbcd = 0
  end subroutine allocate_scalars

  !> @brief Allocate arrays specific to the lake energy transport (LKE)
  !! package.
  !<
  subroutine lke_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweLkeType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- Time series
    call mem_allocate(this%temprain, this%ncv, 'TEMPRAIN', this%memoryPath)
    call mem_allocate(this%tempevap, this%ncv, 'TEMPEVAP', this%memoryPath)
    call mem_allocate(this%temproff, this%ncv, 'TEMPROFF', this%memoryPath)
    call mem_allocate(this%tempiflw, this%ncv, 'TEMPIFLW', this%memoryPath)
    !
    ! -- Call standard TspAptType allocate arrays
    call this%TspAptType%apt_allocate_arrays()
    !
    ! -- Initialize
    do n = 1, this%ncv
      this%temprain(n) = DZERO
      this%tempevap(n) = DZERO
      this%temproff(n) = DZERO
      this%tempiflw(n) = DZERO
    end do
    !
  end subroutine lke_allocate_arrays

  !> @brief Deallocate memory
  !<
  subroutine lke_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GweLkeType) :: this
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%idxbudrain)
    call mem_deallocate(this%idxbudevap)
    call mem_deallocate(this%idxbudroff)
    call mem_deallocate(this%idxbudiflw)
    call mem_deallocate(this%idxbudwdrl)
    call mem_deallocate(this%idxbudoutf)
    call mem_deallocate(this%idxbudlbcd)
    !
    ! -- Deallocate time series
    call mem_deallocate(this%temprain)
    call mem_deallocate(this%tempevap)
    call mem_deallocate(this%temproff)
    call mem_deallocate(this%tempiflw)
    !
    ! -- Deallocate arrays
    call mem_deallocate(this%ktf)
    call mem_deallocate(this%rfeatthk)
    !
    ! -- Deallocate scalars in TspAptType
    call this%TspAptType%bnd_da()
  end subroutine lke_da

  !> @brief Rain term
  !<
  subroutine lke_rain_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    !
    n1 = this%flowbudptr%budterm(this%idxbudrain)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudrain)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudrain)%flow(ientry)
    ctmp = this%temprain(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = -rrate
    if (present(hcofval)) hcofval = DZERO
  end subroutine lke_rain_term

  !> @brief Evaporative term
  !<
  subroutine lke_evap_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: heatlat
    !
    n1 = this%flowbudptr%budterm(this%idxbudevap)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudevap)%id2(ientry)
    ! -- Note that qbnd is negative for evap
    qbnd = this%flowbudptr%budterm(this%idxbudevap)%flow(ientry)
    heatlat = this%gwecommon%gwerhow * this%gwecommon%gwelatheatvap
    if (present(rrate)) rrate = qbnd * heatlat
    if (present(rhsval)) rhsval = -rrate
    if (present(hcofval)) hcofval = DZERO
  end subroutine lke_evap_term

  !> @brief Runoff term
  !<
  subroutine lke_roff_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    !
    n1 = this%flowbudptr%budterm(this%idxbudroff)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudroff)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudroff)%flow(ientry)
    ctmp = this%temproff(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = -rrate
    if (present(hcofval)) hcofval = DZERO
  end subroutine lke_roff_term

  !> @brief Inflow Term
  !!
  !! Accounts for energy flowing into a lake from a connected stream, for
  !! example.
  !<
  subroutine lke_iflw_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    !
    n1 = this%flowbudptr%budterm(this%idxbudiflw)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudiflw)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudiflw)%flow(ientry)
    ctmp = this%tempiflw(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = -rrate
    if (present(hcofval)) hcofval = DZERO
  end subroutine lke_iflw_term

  !> @brief Specified withdrawal term
  !!
  !! Accounts for energy associated with a withdrawal of water from a lake
  !! or group of lakes.
  !<
  subroutine lke_wdrl_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    !
    n1 = this%flowbudptr%budterm(this%idxbudwdrl)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudwdrl)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudwdrl)%flow(ientry)
    ctmp = this%xnewpak(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine lke_wdrl_term

  !> @brief Outflow term
  !!
  !! Accounts for the energy leaving a lake, for example, energy exiting a
  !! lake via a flow into a draining stream channel.
  !<
  subroutine lke_outf_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweLkeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    !
    n1 = this%flowbudptr%budterm(this%idxbudoutf)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudoutf)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudoutf)%flow(ientry)
    ctmp = this%xnewpak(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine lke_outf_term

  !> @brief Defined observation types
  !!
  !! Store the observation type supported by the APT package and override
  !! BndType%bnd_df_obs
  !<
  subroutine lke_df_obs(this)
    ! -- dummy
    class(GweLkeType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for temperature observation type.
    call this%obs%StoreObsType('temperature', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for flow between features, such as lake to lake.
    call this%obs%StoreObsType('flow-ja-face', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID12
    !
    ! -- Store obs type and assign procedure pointer
    !    for from-mvr observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for to-mvr observation type.
    call this%obs%StoreObsType('to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for storage observation type.
    call this%obs%StoreObsType('storage', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for constant observation type.
    call this%obs%StoreObsType('constant', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type: lke
    call this%obs%StoreObsType('lke', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for rainfall observation type.
    call this%obs%StoreObsType('rainfall', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for evaporation observation type.
    call this%obs%StoreObsType('evaporation', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for runoff observation type.
    call this%obs%StoreObsType('runoff', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for inflow observation type.
    call this%obs%StoreObsType('ext-inflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for withdrawal observation type.
    call this%obs%StoreObsType('withdrawal', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for ext-outflow observation type.
    call this%obs%StoreObsType('ext-outflow', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
  end subroutine lke_df_obs

  !> @brief Process package specific obs
  !!
  !! Method to process specific observations for this package.
  !<
  subroutine lke_rp_obs(this, obsrv, found)
    ! -- dummy
    class(GweLkeType), intent(inout) :: this !< package class
    type(ObserveType), intent(inout) :: obsrv !< observation object
    logical, intent(inout) :: found !< indicate whether observation was found
    !
    found = .true.
    select case (obsrv%ObsTypeId)
    case ('RAINFALL')
      call this%rp_obs_byfeature(obsrv)
    case ('EVAPORATION')
      call this%rp_obs_byfeature(obsrv)
    case ('RUNOFF')
      call this%rp_obs_byfeature(obsrv)
    case ('EXT-INFLOW')
      call this%rp_obs_byfeature(obsrv)
    case ('WITHDRAWAL')
      call this%rp_obs_byfeature(obsrv)
    case ('EXT-OUTFLOW')
      call this%rp_obs_byfeature(obsrv)
    case ('TO-MVR')
      call this%rp_obs_budterm(obsrv, &
                               this%flowbudptr%budterm(this%idxbudtmvr))
    case default
      found = .false.
    end select
  end subroutine lke_rp_obs

  !> @brief Calculate observation value and pass it back to APT
  !<
  subroutine lke_bd_obs(this, obstypeid, jj, v, found)
    ! -- dummy
    class(GweLkeType), intent(inout) :: this
    character(len=*), intent(in) :: obstypeid
    real(DP), intent(inout) :: v
    integer(I4B), intent(in) :: jj
    logical, intent(inout) :: found
    ! -- local
    integer(I4B) :: n1, n2
    !
    found = .true.
    select case (obstypeid)
    case ('RAINFALL')
      if (this%iboundpak(jj) /= 0) then
        call this%lke_rain_term(jj, n1, n2, v)
      end if
    case ('EVAPORATION')
      if (this%iboundpak(jj) /= 0) then
        call this%lke_evap_term(jj, n1, n2, v)
      end if
    case ('RUNOFF')
      if (this%iboundpak(jj) /= 0) then
        call this%lke_roff_term(jj, n1, n2, v)
      end if
    case ('EXT-INFLOW')
      if (this%iboundpak(jj) /= 0) then
        call this%lke_iflw_term(jj, n1, n2, v)
      end if
    case ('WITHDRAWAL')
      if (this%iboundpak(jj) /= 0) then
        call this%lke_wdrl_term(jj, n1, n2, v)
      end if
    case ('EXT-OUTFLOW')
      if (this%iboundpak(jj) /= 0) then
        call this%lke_outf_term(jj, n1, n2, v)
      end if
    case default
      found = .false.
    end select
  end subroutine lke_bd_obs

  !> @brief Sets the stress period attributes for keyword use.
  !<
  subroutine lke_set_stressperiod(this, itemno, keyword, found)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(GweLkeType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    character(len=*), intent(in) :: keyword
    logical, intent(inout) :: found
    ! -- local
    character(len=LINELENGTH) :: text
    integer(I4B) :: ierr
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! RAINFALL <rainfall>
    ! EVAPORATION <evaporation>
    ! RUNOFF <runoff>
    ! EXT-INFLOW <inflow>
    ! WITHDRAWAL <withdrawal>
    !
    found = .true.
    select case (keyword)
    case ('RAINFALL')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%temprain(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RAINFALL')
    case ('EVAPORATION')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%tempevap(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'EVAPORATION')
    case ('RUNOFF')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%temproff(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RUNOFF')
    case ('EXT-INFLOW')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%tempiflw(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'EXT-INFLOW')
    case default
      !
      ! -- Keyword not recognized so return to caller with found = .false.
      found = .false.
    end select
    !
999 continue
  end subroutine lke_set_stressperiod

  !> @brief Read feature information for this advanced package
  !<
  subroutine lke_read_cvs(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(GweLkeType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LENBOUNDNAME) :: bndName, bndNameTemp
    character(len=9) :: cno
    character(len=50), dimension(:), allocatable :: caux
    integer(I4B) :: ierr
    logical :: isfound, endOfBlock
    integer(I4B) :: n
    integer(I4B) :: ii, jj
    integer(I4B) :: iaux
    integer(I4B) :: itmp
    integer(I4B) :: nlak
    integer(I4B) :: nconn
    integer(I4B), dimension(:), pointer, contiguous :: nboundchk
    real(DP), pointer :: bndElem => null()
    !
    ! -- initialize itmp
    itmp = 0
    !
    ! -- allocate apt data
    call mem_allocate(this%strt, this%ncv, 'STRT', this%memoryPath)
    call mem_allocate(this%ktf, this%ncv, 'KTF', this%memoryPath)
    call mem_allocate(this%rfeatthk, this%ncv, 'RFEATTHK', this%memoryPath)
    call mem_allocate(this%lauxvar, this%naux, this%ncv, 'LAUXVAR', &
                      this%memoryPath)
    !
    ! -- lake boundary and concentrations
    if (this%imatrows == 0) then
      call mem_allocate(this%iboundpak, this%ncv, 'IBOUND', this%memoryPath)
      call mem_allocate(this%xnewpak, this%ncv, 'XNEWPAK', this%memoryPath)
    end if
    call mem_allocate(this%xoldpak, this%ncv, 'XOLDPAK', this%memoryPath)
    !
    ! -- allocate character storage not managed by the memory manager
    allocate (this%featname(this%ncv)) ! ditch after boundnames allocated??
    !allocate(this%status(this%ncv))
    !
    do n = 1, this%ncv
      this%strt(n) = DEP20
      this%ktf(n) = DZERO
      this%rfeatthk(n) = DZERO
      this%lauxvar(:, n) = DZERO
      this%xoldpak(n) = DEP20
      if (this%imatrows == 0) then
        this%iboundpak(n) = 1
        this%xnewpak(n) = DEP20
      end if
    end do
    !
    ! -- allocate local storage for aux variables
    if (this%naux > 0) then
      allocate (caux(this%naux))
    end if
    !
    ! -- allocate and initialize temporary variables
    allocate (nboundchk(this%ncv))
    do n = 1, this%ncv
      nboundchk(n) = 0
    end do
    !
    ! -- get packagedata block
    call this%parser%GetBlock('PACKAGEDATA', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse locations block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' PACKAGEDATA'
      nlak = 0
      nconn = 0
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        n = this%parser%GetInteger()

        if (n < 1 .or. n > this%ncv) then
          write (errmsg, '(a,1x,i6)') &
            'Itemno must be > 0 and <= ', this%ncv
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- increment nboundchk
        nboundchk(n) = nboundchk(n) + 1
        !
        ! -- strt
        this%strt(n) = this%parser%GetDouble()
        !
        ! -- read additional thermal conductivity terms
        this%ktf(n) = this%parser%GetDouble()
        this%rfeatthk(n) = this%parser%GetDouble()
        if (this%rfeatthk(n) <= DZERO) then
          write (errmsg, '(4x,a)') &
          '****ERROR. Specified thickness used for thermal &
          &conduction MUST BE > 0 else divide by zero error occurs'
          call store_error(errmsg)
          cycle
        end if
        !
        ! -- get aux data
        do iaux = 1, this%naux
          call this%parser%GetString(caux(iaux))
        end do

        ! -- set default bndName
        write (cno, '(i9.9)') n
        bndName = 'Feature'//cno

        ! -- featname
        if (this%inamedbound /= 0) then
          call this%parser%GetStringCaps(bndNameTemp)
          if (bndNameTemp /= '') then
            bndName = bndNameTemp
          end if
        end if
        this%featname(n) = bndName

        ! -- fill time series aware data
        ! -- fill aux data
        do jj = 1, this%naux
          text = caux(jj)
          ii = n
          bndElem => this%lauxvar(jj, ii)
          call read_value_or_time_series_adv(text, ii, jj, bndElem, &
                                             this%packName, 'AUX', &
                                             this%tsManager, this%iprpak, &
                                             this%auxname(jj))
        end do
        !
        nlak = nlak + 1
      end do
      !
      ! -- check for duplicate or missing lakes
      do n = 1, this%ncv
        if (nboundchk(n) == 0) then
          write (errmsg, '(a,1x,i0)') 'No data specified for feature', n
          call store_error(errmsg)
        else if (nboundchk(n) > 1) then
          write (errmsg, '(a,1x,i0,1x,a,1x,i0,1x,a)') &
            'Data for feature', n, 'specified', nboundchk(n), 'times'
          call store_error(errmsg)
        end if
      end do
      !
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' PACKAGEDATA'
    else
      call store_error('Required packagedata block not found.')
    end if
    !
    ! -- terminate if any errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- deallocate local storage for aux variables
    if (this%naux > 0) then
      deallocate (caux)
    end if
    !
    ! -- deallocate local storage for nboundchk
    deallocate (nboundchk)
  end subroutine lke_read_cvs

end module GweLkeModule
