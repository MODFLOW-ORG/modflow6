! -- Multi-Aquifer Well Energy Transport Module
! -- todo: save the mwe temperature into the mwe aux variable?
! -- todo: calculate the maw DENSE aux variable using temperature?
!
! MAW flows (flowbudptr)     index var    MWE term              Transport Type
!---------------------------------------------------------------------------------

! -- terms from MAW that will be handled by parent APT Package
! FLOW-JA-FACE              idxbudfjf     FLOW-JA-FACE          cv2cv  (note that this doesn't exist for MAW)
! GWF (aux FLOW-AREA)       idxbudgwf     GWF                   cv2gwf
! STORAGE (aux VOLUME)      idxbudsto     none                  used for cv volumes
! FROM-MVR                  idxbudfmvr    FROM-MVR              q * text = this%qfrommvr(:)
! TO-MVR                    idxbudtmvr    TO-MVR                q * tfeat

! -- MAW terms
! RATE                      idxbudrate    RATE                  q < 0: q * twell, else q * tuser
! FW-RATE                   idxbudfwrt    FW-RATE               q * twell
! RATE-TO-MVR               idxbudrtmv    RATE-TO-MVR           q * twell
! FW-RATE-TO-MVR            idxbudfrtm    FW-RATE-TO-MVR        q * twell
! WELL-AQUIFER CONDUCTION   idxbudmwcd    MW-CONDUCTION         K_t_f * WetArea / thickness

! -- terms from MAW that should be skipped
! CONSTANT-TO-MVR           ?             CONSTANT-TO-MVR       q * twell

! -- terms from a flow file that should be skipped
! CONSTANT                  none          none                  none
! AUXILIARY                 none          none                  none

! -- terms that are written to the energy transport budget file
! none                      none          STORAGE (aux ENER)    dE/dt
! none                      none          AUXILIARY             none
! none                      none          CONSTANT              accumulate
!
!
module GweMweModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LINELENGTH, LENBOUNDNAME, DEP20
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors
  use BndModule, only: BndType, GetBndFromList
  use TspFmiModule, only: TspFmiType
  use MawModule, only: MawType
  use ObserveModule, only: ObserveType
  use TspAptModule, only: TspAptType, apt_process_obsID, &
                          apt_process_obsID12
  use GweInputDataModule, only: GweInputDataType
  use MatrixBaseModule

  implicit none

  public mwe_create

  character(len=*), parameter :: ftype = 'MWE'
  character(len=*), parameter :: flowtype = 'MAW'
  character(len=16) :: text = '             MWE'

  type, extends(TspAptType) :: GweMweType

    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in mst

    integer(I4B), pointer :: idxbudrate => null() ! index of well rate terms in flowbudptr
    integer(I4B), pointer :: idxbudfwrt => null() ! index of flowing well rate terms in flowbudptr
    integer(I4B), pointer :: idxbudrtmv => null() ! index of rate to mover terms in flowbudptr
    integer(I4B), pointer :: idxbudfrtm => null() ! index of flowing well rate to mover terms in flowbudptr
    integer(I4B), pointer :: idxbudmwcd => null() ! index of well bore conduction terms in flowbudptr
    real(DP), dimension(:), pointer, contiguous :: ktf => null() !< thermal conductivity between the sfe and groundwater cell
    real(DP), dimension(:), pointer, contiguous :: rfeatthk => null() !< thickness of streambed material through which thermal conduction occurs
    real(DP), dimension(:), pointer, contiguous :: temprate => null() ! well rate temperature

  contains

    procedure :: bnd_da => mwe_da
    procedure :: allocate_scalars
    procedure :: apt_allocate_arrays => mwe_allocate_arrays
    procedure :: find_apt_package => find_mwe_package
    procedure :: pak_fc_expanded => mwe_fc_expanded
    procedure :: pak_solve => mwe_solve
    procedure :: pak_get_nbudterms => mwe_get_nbudterms
    procedure :: pak_setup_budobj => mwe_setup_budobj
    procedure :: pak_fill_budobj => mwe_fill_budobj
    procedure :: mwe_rate_term
    procedure :: mwe_fwrt_term
    procedure :: mwe_rtmv_term
    procedure :: mwe_frtm_term
    procedure :: pak_df_obs => mwe_df_obs
    procedure :: pak_rp_obs => mwe_rp_obs
    procedure :: pak_bd_obs => mwe_bd_obs
    procedure :: pak_set_stressperiod => mwe_set_stressperiod
    procedure :: apt_read_cvs => mwe_read_cvs

  end type GweMweType

contains

  !> Create new MWE package
  !<
  subroutine mwe_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(GweMweType), pointer :: mweobj
    !
    ! -- Allocate the object and assign values to object variables
    allocate (mweobj)
    packobj => mweobj
    !
    ! -- Create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- Allocate scalars
    call mweobj%allocate_scalars()
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
    mweobj%fmi => fmi
    !
    ! -- Store pointer to governing equation scale factor
    mweobj%eqnsclfac => eqnsclfac
    !
    ! -- Store pointer to shared data module for accessing cpw, rhow
    !    for the budget calculations, and for accessing the latent heat of
    !    vaporization for evaporative cooling.
    mweobj%gwecommon => gwecommon
    !
    ! -- Set labels that will be used in generalized APT class
    mweobj%depvartype = dvt
    mweobj%depvarunit = dvu
    mweobj%depvarunitabbrev = dvua
  end subroutine mwe_create

  !> @brief Find corresponding mwe package
  !<
  subroutine find_mwe_package(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweMweType) :: this
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
            type is (MawType)
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
        this%idxbudssm(ip) = 0
      case ('STORAGE')
        this%idxbudsto = ip
        this%idxbudssm(ip) = 0
      case ('RATE')
        this%idxbudrate = ip
        this%idxbudssm(ip) = 0
      case ('FW-RATE')
        this%idxbudfwrt = ip
        this%idxbudssm(ip) = 0
      case ('RATE-TO-MVR')
        this%idxbudrtmv = ip
        this%idxbudssm(ip) = 0
      case ('FW-RATE-TO-MVR')
        this%idxbudfrtm = ip
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
    !
    ! -- Streambed conduction term
    this%idxbudmwcd = this%idxbudgwf
  end subroutine find_mwe_package

  !> @brief Add matrix terms related to MWE
  !!
  !! This routine is called from TspAptType%apt_fc_expanded() in
  !! order to add matrix terms specifically for MWE
  !<
  subroutine mwe_fc_expanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GweMweType) :: this
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
    real(DP) :: ctherm ! kluge?
    real(DP) :: wa !< wetted area
    real(DP) :: ktf !< thermal conductivity of streambed material
    real(DP) :: s !< thickness of conductive wellbore material
    !
    ! -- Add puping rate contribution
    if (this%idxbudrate /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrate)%nlist
        call this%mwe_rate_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add flowing well rate contribution
    if (this%idxbudfwrt /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfwrt)%nlist
        call this%mwe_fwrt_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add rate to mover contribution
    if (this%idxbudrtmv /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrtmv)%nlist
        call this%mwe_rtmv_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add puping rate contribution
    if (this%idxbudfrtm /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfrtm)%nlist
        call this%mwe_frtm_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add wellbore conduction contribution
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      !
      ! -- Set n to feature number and process if active features
      n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      if (this%iboundpak(n) /= 0) then
        !
        ! -- Set acoef and rhs to negative so they are relative to mwe and not gwe
        auxpos = this%flowbudptr%budterm(this%idxbudgwf)%naux
        wa = this%flowbudptr%budterm(this%idxbudgwf)%auxvar(auxpos, j)
        ktf = this%ktf(n)
        s = this%rfeatthk(n)
        ctherm = ktf * wa / s
        !
        ! -- Add to mwe row
        iposd = this%idxdglo(j)
        iposoffd = this%idxoffdglo(j)
        call matrix_sln%add_value_pos(iposd, -ctherm) ! kluge note: make sure the signs on ctherm are correct here and below
        call matrix_sln%add_value_pos(iposoffd, ctherm)
        !
        ! -- Add to gwe row for mwe connection
        ipossymd = this%idxsymdglo(j)
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymd, -ctherm)
        call matrix_sln%add_value_pos(ipossymoffd, ctherm)
      end if
    end do
  end subroutine mwe_fc_expanded

  !> @brief Add terms specific to multi-aquifer wells to the explicit multi-
  !! aquifer well energy transport solve
  !<
  subroutine mwe_solve(this)
    ! -- dummy
    class(GweMweType) :: this
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n1, n2
    real(DP) :: rrate
    !
    ! -- Add well pumping contribution
    if (this%idxbudrate /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrate)%nlist
        call this%mwe_rate_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add flowing well rate contribution
    if (this%idxbudfwrt /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfwrt)%nlist
        call this%mwe_fwrt_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add well pumping rate to mover contribution
    if (this%idxbudrtmv /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrtmv)%nlist
        call this%mwe_rtmv_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add flowing well rate to mover contribution
    if (this%idxbudfrtm /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfrtm)%nlist
        call this%mwe_frtm_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
  end subroutine mwe_solve

  !> @brief Function to return the number of budget terms just for this package
  !!
  !! This overrides a function in the parent class.
  !<
  function mwe_get_nbudterms(this) result(nbudterms)
    ! -- dummy
    class(GweMweType) :: this
    ! -- return
    integer(I4B) :: nbudterms
    !
    ! -- Number of potential budget terms is 5
    nbudterms = 1 ! RATE
    if (this%idxbudfwrt /= 0) nbudterms = nbudterms + 1
    if (this%idxbudrtmv /= 0) nbudterms = nbudterms + 1
    if (this%idxbudfrtm /= 0) nbudterms = nbudterms + 1
    if (this%idxbudmwcd /= 0) nbudterms = nbudterms + 1
  end function mwe_get_nbudterms

  !> @brief Set up the budget object that stores all the mwe flows
  !<
  subroutine mwe_setup_budobj(this, idx)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(GweMweType) :: this
    integer(I4B), intent(inout) :: idx
    ! -- local
    integer(I4B) :: n, n1, n2
    integer(I4B) :: maxlist, naux
    real(DP) :: q
    character(len=LENBUDTXT) :: text
    !
    ! -- User-specified rate
    text = '            RATE'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudrate)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Flowing well rate
    if (this%idxbudfwrt /= 0) then
      text = '         FW-RATE'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudfwrt)%maxlist
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
    end if
    !
    ! -- User-specified flow rate to mover
    if (this%idxbudrtmv /= 0) then
      text = '     RATE-TO-MVR'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudrtmv)%maxlist
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
    end if
    !
    ! -- Fflowing well rate to mover
    if (this%idxbudfrtm /= 0) then
      text = '  FW-RATE-TO-MVR'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudfrtm)%maxlist
      naux = 0
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux)
    end if
    !
    ! -- Conduction through wellbore (and/or filter pack)
    text = '   WELLBORE-COND'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudmwcd)%maxlist
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
  end subroutine mwe_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  subroutine mwe_fill_budobj(this, idx, x, flowja, ccratin, ccratout)
    ! -- dummy
    class(GweMweType) :: this
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
    real(DP) :: ctherm
    real(DP) :: wa !< wetted area
    real(DP) :: ktf !< thermal conductivity of streambed material
    real(DP) :: s !< thickness of conductive streambed materia
    !
    ! -- Rate
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudrate)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%mwe_rate_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- FW-Rate
    if (this%idxbudfwrt /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudfwrt)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%mwe_fwrt_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- Rate-To-MVR
    if (this%idxbudrtmv /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudrtmv)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%mwe_rtmv_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- FW-Rate-To-MVR
    if (this%idxbudfrtm /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudfrtm)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%mwe_frtm_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- Wellbore-Cond
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do j = 1, this%flowbudptr%budterm(this%idxbudmwcd)%nlist
      q = DZERO
      n1 = this%flowbudptr%budterm(this%idxbudmwcd)%id1(j)
      if (this%iboundpak(n1) /= 0) then
        igwfnode = this%flowbudptr%budterm(this%idxbudmwcd)%id2(j)
        auxpos = this%flowbudptr%budterm(this%idxbudgwf)%naux
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
        this%simvals(j) = this%simvals(j) - q
        idiag = this%dis%con%ia(igwfnode)
        flowja(idiag) = flowja(idiag) - q
      end if
    end do
  end subroutine mwe_fill_budobj

  !> @brief Allocate scalars specific to the multi-aquifer well energy
  !! transport (MWE) package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweMweType) :: this
    !
    ! -- Allocate scalars in TspAptType
    call this%TspAptType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%idxbudrate, 'IDXBUDRATE', this%memoryPath)
    call mem_allocate(this%idxbudfwrt, 'IDXBUDFWRT', this%memoryPath)
    call mem_allocate(this%idxbudrtmv, 'IDXBUDRTMV', this%memoryPath)
    call mem_allocate(this%idxbudfrtm, 'IDXBUDFRTM', this%memoryPath)
    call mem_allocate(this%idxbudmwcd, 'IDXBUDMWCD', this%memoryPath)
    !
    ! -- Initialize
    this%idxbudrate = 0
    this%idxbudfwrt = 0
    this%idxbudrtmv = 0
    this%idxbudfrtm = 0
    this%idxbudmwcd = 0
  end subroutine allocate_scalars

  !> @brief Allocate arrays specific to the streamflow mass transport (SFT)
  !! package
  !<
  subroutine mwe_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweMweType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- Time series
    call mem_allocate(this%temprate, this%ncv, 'TEMPRATE', this%memoryPath)
    !
    ! -- Call standard TspAptType allocate arrays
    call this%TspAptType%apt_allocate_arrays()
    !
    ! -- Initialize
    do n = 1, this%ncv
      this%temprate(n) = DZERO
    end do
  end subroutine mwe_allocate_arrays

  !> @brief Deallocate memory associated with MWE package
  !<
  subroutine mwe_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GweMweType) :: this
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%idxbudrate)
    call mem_deallocate(this%idxbudfwrt)
    call mem_deallocate(this%idxbudrtmv)
    call mem_deallocate(this%idxbudfrtm)
    call mem_deallocate(this%idxbudmwcd)
    !
    ! -- Deallocate time series
    call mem_deallocate(this%temprate)
    !
    ! -- Deallocate arrays
    call mem_deallocate(this%ktf)
    call mem_deallocate(this%rfeatthk)
    !
    ! -- Deallocate scalars in TspAptType
    call this%TspAptType%bnd_da()
  end subroutine mwe_da

  !> @brief Thermal transport matrix term(s) associated with a user-specified
  !! flow rate (mwe_rate_term)
  !<
  subroutine mwe_rate_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(GweMweType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    real(DP) :: h, r
    !
    n1 = this%flowbudptr%budterm(this%idxbudrate)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudrate)%id2(ientry)
    ! -- Note that qbnd is negative for an extracting well
    qbnd = this%flowbudptr%budterm(this%idxbudrate)%flow(ientry)
    if (qbnd < DZERO) then
      ctmp = this%xnewpak(n1)
      h = qbnd
      r = DZERO
    else
      ctmp = this%temprate(n1)
      h = DZERO
      r = -qbnd * ctmp
    end if
    if (present(rrate)) rrate = qbnd * ctmp * this%eqnsclfac
    if (present(rhsval)) rhsval = r * this%eqnsclfac
    if (present(hcofval)) hcofval = h * this%eqnsclfac
  end subroutine mwe_rate_term

  !> @brief Thermal transport matrix term(s) associated with a flowing-
  !! well rate term associated with pumping (or injection)
  !<
  subroutine mwe_fwrt_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(GweMweType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudfwrt)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudfwrt)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudfwrt)%flow(ientry)
    ctmp = this%xnewpak(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine mwe_fwrt_term

  !> @brief Thermal transport matrix term(s) associated with pumped-water-
  !! to-mover term (mwe_rtmv_term)
  !<
  subroutine mwe_rtmv_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(GweMweType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudrtmv)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudrtmv)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudrtmv)%flow(ientry)
    ctmp = this%xnewpak(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine mwe_rtmv_term

  !> @brief Thermal transport matrix term(s) associated with the flowing-
  !! well-rate-to-mover term (mwe_frtm_term)
  !<
  subroutine mwe_frtm_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(GweMweType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudfrtm)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudfrtm)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudfrtm)%flow(ientry)
    ctmp = this%xnewpak(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine mwe_frtm_term

  !> @brief Observations
  !!
  !! Store the observation type supported by the APT package and override
  !! BndType%bnd_df_obs
  !<
  subroutine mwe_df_obs(this)
    ! -- dummy
    class(GweMweType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for temperature observation type.
    call this%obs%StoreObsType('temperature', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Flow-ja-face not supported for MWE
    !call this%obs%StoreObsType('flow-ja-face', .true., indx)
    !this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for from-mvr observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- To-mvr not supported for mwe
    !call this%obs%StoreObsType('to-mvr', .true., indx)
    !this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
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
    !    for observation type: mwe
    call this%obs%StoreObsType('mwe', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID12
    !
    ! -- Store obs type and assign procedure pointer
    !    for rate observation type.
    call this%obs%StoreObsType('rate', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('fw-rate', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('rate-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('fw-rate-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
  end subroutine mwe_df_obs

  !> @brief Process package specific obs
  !!
  !! Method to process specific observations for this package.
  !<
  subroutine mwe_rp_obs(this, obsrv, found)
    ! -- dummy
    class(GweMweType), intent(inout) :: this !< package class
    type(ObserveType), intent(inout) :: obsrv !< observation object
    logical, intent(inout) :: found !< indicate whether observation was found
    !
    found = .true.
    select case (obsrv%ObsTypeId)
    case ('RATE')
      call this%rp_obs_byfeature(obsrv)
    case ('FW-RATE')
      call this%rp_obs_byfeature(obsrv)
    case ('RATE-TO-MVR')
      call this%rp_obs_byfeature(obsrv)
    case ('FW-RATE-TO-MVR')
      call this%rp_obs_byfeature(obsrv)
    case default
      found = .false.
    end select
  end subroutine mwe_rp_obs

  !> @brief Calculate observation value and pass it back to APT
  !<
  subroutine mwe_bd_obs(this, obstypeid, jj, v, found)
    ! -- dummy
    class(GweMweType), intent(inout) :: this
    character(len=*), intent(in) :: obstypeid
    real(DP), intent(inout) :: v
    integer(I4B), intent(in) :: jj
    logical, intent(inout) :: found
    ! -- local
    integer(I4B) :: n1, n2
    !
    found = .true.
    select case (obstypeid)
    case ('RATE')
      if (this%iboundpak(jj) /= 0) then
        call this%mwe_rate_term(jj, n1, n2, v)
      end if
    case ('FW-RATE')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudfwrt > 0) then
        call this%mwe_fwrt_term(jj, n1, n2, v)
      end if
    case ('RATE-TO-MVR')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudrtmv > 0) then
        call this%mwe_rtmv_term(jj, n1, n2, v)
      end if
    case ('FW-RATE-TO-MVR')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudfrtm > 0) then
        call this%mwe_frtm_term(jj, n1, n2, v)
      end if
    case default
      found = .false.
    end select
  end subroutine mwe_bd_obs

  !> @brief Sets the stress period attributes for keyword use.
  !<
  subroutine mwe_set_stressperiod(this, itemno, keyword, found)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(GweMweType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    character(len=*), intent(in) :: keyword
    logical, intent(inout) :: found
    ! -- local
    character(len=LINELENGTH) :: text
    integer(I4B) :: ierr
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! RATE <rate>
    !
    found = .true.
    select case (keyword)
    case ('RATE')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1
      bndElem => this%temprate(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RATE')
    case default
      !
      ! -- Keyword not recognized so return to caller with found = .false.
      found = .false.
    end select
    !
999 continue
  end subroutine mwe_set_stressperiod

  !> @brief Read feature information for this advanced package
  !<
  subroutine mwe_read_cvs(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(GweMweType), intent(inout) :: this
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
  end subroutine mwe_read_cvs

end module GweMweModule
