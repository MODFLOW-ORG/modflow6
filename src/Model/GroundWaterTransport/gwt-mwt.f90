! -- Multi-Aquifer Well Transport Module
! -- todo: what to do about reactions in maw?  Decay?
! -- todo: save the mwt concentration into the mwt aux variable?
! -- todo: calculate the maw DENSE aux variable using concentration?
!
! MAW flows (flowbudptr)     index var    MWT term              Transport Type
!---------------------------------------------------------------------------------

! -- terms from MAW that will be handled by parent APT Package
! FLOW-JA-FACE              idxbudfjf     FLOW-JA-FACE          cv2cv  (note that this doesn't exist for MAW)
! GWF (aux FLOW-AREA)       idxbudgwf     GWF                   cv2gwf
! STORAGE (aux VOLUME)      idxbudsto     none                  used for cv volumes
! FROM-MVR                  idxbudfmvr    FROM-MVR              q * cext = this%qfrommvr(:)
! TO-MVR                    idxbudtmvr    TO-MVR                q * cfeat

! -- MAW terms
! RATE                      idxbudrate    RATE                  q < 0: q * cwell, else q * cuser
! FW-RATE                   idxbudfwrt    FW-RATE               q * cwell
! RATE-TO-MVR               idxbudrtmv    RATE-TO-MVR           q * cwell
! FW-RATE-TO-MVR            idxbudfrtm    FW-RATE-TO-MVR        q * cwell

! -- terms from MAW that should be skipped
! CONSTANT-TO-MVR           ?             CONSTANT-TO-MVR       q * cwell

! -- terms from a flow file that should be skipped
! CONSTANT                  none          none                  none
! AUXILIARY                 none          none                  none

! -- terms that are written to the transport budget file
! none                      none          STORAGE (aux MASS)    dM/dt
! none                      none          AUXILIARY             none
! none                      none          CONSTANT              accumulate
!
!
module GwtMwtModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, LINELENGTH
  use SimModule, only: store_error
  use BndModule, only: BndType, GetBndFromList
  use TspFmiModule, only: TspFmiType
  use MawModule, only: MawType
  use ObserveModule, only: ObserveType
  use TspAptModule, only: TspAptType, apt_process_obsID, &
                          apt_process_obsID12
  use MatrixBaseModule

  implicit none

  public mwt_create

  character(len=*), parameter :: ftype = 'MWT'
  character(len=*), parameter :: flowtype = 'MAW'
  character(len=16) :: text = '             MWT'

  type, extends(TspAptType) :: GwtMwtType

    integer(I4B), pointer :: idxbudrate => null() ! index of well rate terms in flowbudptr
    integer(I4B), pointer :: idxbudfwrt => null() ! index of flowing well rate terms in flowbudptr
    integer(I4B), pointer :: idxbudrtmv => null() ! index of rate to mover terms in flowbudptr
    integer(I4B), pointer :: idxbudfrtm => null() ! index of flowing well rate to mover terms in flowbudptr
    real(DP), dimension(:), pointer, contiguous :: concrate => null() ! well rate concentration

  contains

    procedure :: bnd_da => mwt_da
    procedure :: allocate_scalars
    procedure :: apt_allocate_arrays => mwt_allocate_arrays
    procedure :: find_apt_package => find_mwt_package
    procedure :: pak_fc_expanded => mwt_fc_expanded
    procedure :: pak_solve => mwt_solve
    procedure :: pak_get_nbudterms => mwt_get_nbudterms
    procedure :: pak_setup_budobj => mwt_setup_budobj
    procedure :: pak_fill_budobj => mwt_fill_budobj
    procedure :: mwt_rate_term
    procedure :: mwt_fwrt_term
    procedure :: mwt_rtmv_term
    procedure :: mwt_frtm_term
    procedure :: pak_df_obs => mwt_df_obs
    procedure :: pak_rp_obs => mwt_rp_obs
    procedure :: pak_bd_obs => mwt_bd_obs
    procedure :: pak_set_stressperiod => mwt_set_stressperiod

  end type GwtMwtType

contains

  !> Create new MWT package
  !<
  subroutine mwt_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        fmi, eqnsclfac, dvt, dvu, dvua)
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
    character(len=*), intent(in) :: dvt !< For GWT, set to "CONCENTRATION" in TspAptType
    character(len=*), intent(in) :: dvu !< For GWT, set to "mass" in TspAptType
    character(len=*), intent(in) :: dvua !< For GWT, set to "M" in TspAptType
    ! -- local
    type(GwtMwtType), pointer :: mwtobj
    !
    ! -- allocate the object and assign values to object variables
    allocate (mwtobj)
    packobj => mwtobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- allocate scalars
    call mwtobj%allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ncolbnd = 1
    packobj%iscloc = 1

    ! -- Store pointer to flow model interface.  When the GwfGwt exchange is
    !    created, it sets fmi%bndlist so that the GWT model has access to all
    !    the flow packages
    mwtobj%fmi => fmi
    !
    ! -- Store pointer to governing equation scale factor
    mwtobj%eqnsclfac => eqnsclfac
    !
    ! -- Set labels that will be used in generalized APT class
    mwtobj%depvartype = dvt
    mwtobj%depvarunit = dvu
    mwtobj%depvarunitabbrev = dvua
  end subroutine mwt_create

  !> @brief find corresponding mwt package
  !<
  subroutine find_mwt_package(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtMwtType) :: this
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
            ! -- store BndType pointer to packobj, and then
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
    ! -- error if flow package not found
    if (.not. found) then
      write (errmsg, '(a)') 'Could not find flow package with name '&
                            &//trim(adjustl(this%flowpackagename))//'.'
      call store_error(errmsg)
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- allocate space for idxbudssm, which indicates whether this is a
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
        ! -- set idxbudssm equal to a column index for where the concentrations
        !    are stored in the concbud(nbudssm, ncv) array
        this%idxbudssm(ip) = icount
        icount = icount + 1
      end select
      write (this%iout, '(a, i0, " = ", a,/, a, i0)') &
        '  TERM ', ip, trim(adjustl(this%flowbudptr%budterm(ip)%flowtype)), &
        '   MAX NO. OF ENTRIES = ', this%flowbudptr%budterm(ip)%maxlist
    end do
    write (this%iout, '(a, //)') 'DONE PROCESSING '//ftype//' INFORMATION'
  end subroutine find_mwt_package

  !> @brief Add matrix terms related to MWT
  !!
  !! This routine is called from TspAptType%apt_fc_expanded() in
  !! order to add matrix terms specifically for MWT
  !<
  subroutine mwt_fc_expanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- modules
    ! -- dummy
    class(GwtMwtType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, n1, n2
    integer(I4B) :: iloc
    integer(I4B) :: iposd
    real(DP) :: rrate
    real(DP) :: rhsval
    real(DP) :: hcofval
    !
    ! -- add puping rate contribution
    if (this%idxbudrate /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrate)%nlist
        call this%mwt_rate_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- add flowing well rate contribution
    if (this%idxbudfwrt /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfwrt)%nlist
        call this%mwt_fwrt_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- add rate to mover contribution
    if (this%idxbudrtmv /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrtmv)%nlist
        call this%mwt_rtmv_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- add puping rate contribution
    if (this%idxbudfrtm /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfrtm)%nlist
        call this%mwt_frtm_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
  end subroutine mwt_fc_expanded

  !> @ brief Add terms specific to multi-aquifer wells to the explicit multi-
  !! aquifer well solute transport solve
  !<
  subroutine mwt_solve(this)
    ! -- dummy
    class(GwtMwtType) :: this
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n1, n2
    real(DP) :: rrate
    !
    ! -- add well pumping contribution
    if (this%idxbudrate /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrate)%nlist
        call this%mwt_rate_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- add flowing well rate contribution
    if (this%idxbudfwrt /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfwrt)%nlist
        call this%mwt_fwrt_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- add well pumping rate to mover contribution
    if (this%idxbudrtmv /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrtmv)%nlist
        call this%mwt_rtmv_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- add flowing well rate to mover contribution
    if (this%idxbudfrtm /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfrtm)%nlist
        call this%mwt_frtm_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
  end subroutine mwt_solve

  !> @brief Function to return the number of budget terms just for this package
  !!
  !! This overrides a function in the parent class.
  !<
  function mwt_get_nbudterms(this) result(nbudterms)
    ! -- modules
    ! -- dummy
    class(GwtMwtType) :: this
    ! -- Return
    integer(I4B) :: nbudterms
    ! -- local
    !
    ! -- Number of budget terms is 4
    nbudterms = 1
    if (this%idxbudfwrt /= 0) nbudterms = nbudterms + 1
    if (this%idxbudrtmv /= 0) nbudterms = nbudterms + 1
    if (this%idxbudfrtm /= 0) nbudterms = nbudterms + 1
  end function mwt_get_nbudterms

  !> @brief Set up the budget object that stores all the mwt flows
  !<
  subroutine mwt_setup_budobj(this, idx)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(GwtMwtType) :: this
    integer(I4B), intent(inout) :: idx
    ! -- local
    integer(I4B) :: maxlist, naux
    character(len=LENBUDTXT) :: text
    !
    ! --
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
    ! --
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
    ! --
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
    ! --
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
  end subroutine mwt_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  subroutine mwt_fill_budobj(this, idx, x, flowja, ccratin, ccratout)
    ! -- modules
    ! -- dummy
    class(GwtMwtType) :: this
    integer(I4B), intent(inout) :: idx
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    real(DP), intent(inout) :: ccratin
    real(DP), intent(inout) :: ccratout
    ! -- local
    integer(I4B) :: j, n1, n2
    integer(I4B) :: nlist
    real(DP) :: q
    ! -- formats
    !
    ! -- RATE
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudrate)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%mwt_rate_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- FW-RATE
    if (this%idxbudfwrt /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudfwrt)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%mwt_fwrt_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- RATE-TO-MVR
    if (this%idxbudrtmv /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudrtmv)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%mwt_rtmv_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- FW-RATE-TO-MVR
    if (this%idxbudfrtm /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudfrtm)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%mwt_frtm_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
  end subroutine mwt_fill_budobj

  !> @brief Allocate scalars specific to the streamflow mass transport (SFT)
  !! package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtMwtType) :: this
    ! -- local
    !
    ! -- allocate scalars in TspAptType
    call this%TspAptType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%idxbudrate, 'IDXBUDRATE', this%memoryPath)
    call mem_allocate(this%idxbudfwrt, 'IDXBUDFWRT', this%memoryPath)
    call mem_allocate(this%idxbudrtmv, 'IDXBUDRTMV', this%memoryPath)
    call mem_allocate(this%idxbudfrtm, 'IDXBUDFRTM', this%memoryPath)
    !
    ! -- Initialize
    this%idxbudrate = 0
    this%idxbudfwrt = 0
    this%idxbudrtmv = 0
    this%idxbudfrtm = 0
  end subroutine allocate_scalars

  !> @brief Allocate arrays specific to the streamflow mass transport (SFT)
  !! package.
  !<
  subroutine mwt_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwtMwtType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- time series
    call mem_allocate(this%concrate, this%ncv, 'CONCRATE', this%memoryPath)
    !
    ! -- call standard TspAptType allocate arrays
    call this%TspAptType%apt_allocate_arrays()
    !
    ! -- Initialize
    do n = 1, this%ncv
      this%concrate(n) = DZERO
    end do
    !
  end subroutine mwt_allocate_arrays

  !> @brief Deallocate memory
  !<
  subroutine mwt_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwtMwtType) :: this
    ! -- local
    !
    ! -- deallocate scalars
    call mem_deallocate(this%idxbudrate)
    call mem_deallocate(this%idxbudfwrt)
    call mem_deallocate(this%idxbudrtmv)
    call mem_deallocate(this%idxbudfrtm)
    !
    ! -- deallocate time series
    call mem_deallocate(this%concrate)
    !
    ! -- deallocate scalars in TspAptType
    call this%TspAptType%bnd_da()
  end subroutine mwt_da

  !> @brief Rate term associated with pumping (or injection)
  !<
  subroutine mwt_rate_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GwtMwtType) :: this
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
    ! -- note that qbnd is negative for extracting well
    qbnd = this%flowbudptr%budterm(this%idxbudrate)%flow(ientry)
    if (qbnd < DZERO) then
      ctmp = this%xnewpak(n1)
      h = qbnd
      r = DZERO
    else
      ctmp = this%concrate(n1)
      h = DZERO
      r = -qbnd * ctmp
    end if
    if (present(rrate)) rrate = qbnd * ctmp
    if (present(rhsval)) rhsval = r
    if (present(hcofval)) hcofval = h
  end subroutine mwt_rate_term

  !> @brief Transport matrix term(s) associated with a flowing-
  !! well rate term associated with pumping (or injection)
  !<
  subroutine mwt_fwrt_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GwtMwtType) :: this
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
    if (present(rrate)) rrate = ctmp * qbnd
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd
  end subroutine mwt_fwrt_term

  !> @brief Rate-to-mvr term associated with pumping (or injection)
  !!
  !! Pumped water that is made available to the MVR package for transfer to
  !! another advanced package
  !<
  subroutine mwt_rtmv_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GwtMwtType) :: this
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
    if (present(rrate)) rrate = ctmp * qbnd
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd
  end subroutine mwt_rtmv_term

  !> @brief Flowing well rate-to-mvr term (or injection)
  !!
  !! Pumped water that is made available to the MVR package for transfer to
  !! another advanced package
  !<
  subroutine mwt_frtm_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GwtMwtType) :: this
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
    if (present(rrate)) rrate = ctmp * qbnd
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd
  end subroutine mwt_frtm_term

  !> @brief Observations
  !!
  !! Store the observation type supported by the APT package and override
  !! BndType%bnd_df_obs
  !<
  subroutine mwt_df_obs(this)
    ! -- modules
    ! -- dummy
    class(GwtMwtType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for concentration observation type.
    call this%obs%StoreObsType('concentration', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- flow-ja-face not supported for MWT
    !call this%obs%StoreObsType('flow-ja-face', .true., indx)
    !this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for from-mvr observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- to-mvr not supported for mwt
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
    !    for observation type: mwt
    call this%obs%StoreObsType('mwt', .true., indx)
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
  end subroutine mwt_df_obs

  !> @brief Process package specific obs
  !!
  !! Method to process specific observations for this package.
  !<
  subroutine mwt_rp_obs(this, obsrv, found)
    ! -- dummy
    class(GwtMwtType), intent(inout) :: this !< package class
    type(ObserveType), intent(inout) :: obsrv !< observation object
    logical, intent(inout) :: found !< indicate whether observation was found
    ! -- local
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
  end subroutine mwt_rp_obs

  !> @brief Calculate observation value and pass it back to APT
  !<
  subroutine mwt_bd_obs(this, obstypeid, jj, v, found)
    ! -- dummy
    class(GwtMwtType), intent(inout) :: this
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
        call this%mwt_rate_term(jj, n1, n2, v)
      end if
    case ('FW-RATE')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudfwrt > 0) then
        call this%mwt_fwrt_term(jj, n1, n2, v)
      end if
    case ('RATE-TO-MVR')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudrtmv > 0) then
        call this%mwt_rtmv_term(jj, n1, n2, v)
      end if
    case ('FW-RATE-TO-MVR')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudfrtm > 0) then
        call this%mwt_frtm_term(jj, n1, n2, v)
      end if
    case default
      found = .false.
    end select
  end subroutine mwt_bd_obs

  !> @brief Sets the stress period attributes for keyword use.
  !<
  subroutine mwt_set_stressperiod(this, itemno, keyword, found)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(GwtMwtType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    character(len=*), intent(in) :: keyword
    logical, intent(inout) :: found
    ! -- local
    character(len=LINELENGTH) :: text
    integer(I4B) :: ierr
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    ! -- formats
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
      bndElem => this%concrate(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'RATE')
    case default
      !
      ! -- keyword not recognized so return to caller with found = .false.
      found = .false.
    end select
    !
999 continue
  end subroutine mwt_set_stressperiod

end module GwtMwtModule
