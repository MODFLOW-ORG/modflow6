! -- Unsaturated Zone Flow Energy Transport Module
! -- todo: save the uze temperature into the uze aux variable?
! -- todo: calculate the uzf DENSE aux variable using temperature?
! -- todo: GWD and GWD-TO-MVR do not seem to be included; prob in UZF?
!
! UZF flows (flowbudptr)     index var    UZE term              Transport Type
!---------------------------------------------------------------------------------

! -- terms from UZF that will be handled by parent APT Package
! FLOW-JA-FACE              idxbudfjf     FLOW-JA-FACE          cv2cv
! GWF (aux FLOW-AREA)       idxbudgwf     GWF                   uzf2gwf
! STORAGE (aux VOLUME)      idxbudsto     none                  used for water volumes
! FROM-MVR                  idxbudfmvr    FROM-MVR              q * cext = this%qfrommvr(:)
! AUXILIARY                 none          none                  none
! none                      none          STORAGE (aux MASS)
! none                      none          AUXILIARY             none

! -- terms from UZF that need to be handled here
! INFILTRATION              idxbudinfl    INFILTRATION          q < 0: q * cwell, else q * cuser
! REJ-INF                   idxbudrinf    REJ-INF               q * cuze
! UZET                      idxbuduzet    UZET                  q * cet
! REJ-INF-TO-MVR            idxbudritm    REJ-INF-TO-MVR        q * cinfil?
! THERMAL-EQUIL             idxbudtheq    THERMAL-EQUIL         residual

! -- terms from UZF that should be skipped

module GweUzeModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: DZERO, DONE, LINELENGTH
  use SimModule, only: store_error
  use BndModule, only: BndType, GetBndFromList
  use TspFmiModule, only: TspFmiType
  use UzfModule, only: UzfType
  use ObserveModule, only: ObserveType
  use TspAptModule, only: TspAptType, apt_process_obsID, &
                          apt_process_obsID12
  use GweInputDataModule, only: GweInputDataType
  use MatrixBaseModule

  implicit none

  public uze_create

  character(len=*), parameter :: ftype = 'UZE'
  character(len=*), parameter :: flowtype = 'UZF'
  character(len=16) :: text = '             UZE'

  type, extends(TspAptType) :: GweUzeType

    type(GweInputDataType), pointer :: gwecommon => null() !< pointer to shared gwe data used by multiple packages but set in mst

    integer(I4B), pointer :: idxbudinfl => null() !< index of uzf infiltration terms in flowbudptr
    integer(I4B), pointer :: idxbudrinf => null() !< index of rejected infiltration terms in flowbudptr
    integer(I4B), pointer :: idxbuduzet => null() !< index of unsat et terms in flowbudptr
    integer(I4B), pointer :: idxbudritm => null() !< index of rej infil to mover rate to mover terms in flowbudptr
    integer(I4B), pointer :: idxbudtheq => null() !< index of thermal equilibration terms in flowbudptr

    real(DP), dimension(:), pointer, contiguous :: tempinfl => null() !< infiltration temperature
    real(DP), dimension(:), pointer, contiguous :: tempuzet => null() !< unsat et temperature

  contains

    procedure :: bnd_da => uze_da
    procedure :: allocate_scalars
    procedure :: apt_allocate_arrays => uze_allocate_arrays
    procedure :: find_apt_package => find_uze_package
    procedure :: apt_fc_expanded => uze_fc_expanded
    procedure :: pak_solve => uze_solve
    procedure :: pak_get_nbudterms => uze_get_nbudterms
    procedure :: pak_setup_budobj => uze_setup_budobj
    procedure :: pak_fill_budobj => uze_fill_budobj
    procedure :: uze_infl_term
    procedure :: uze_rinf_term
    procedure :: uze_uzet_term
    procedure :: uze_ritm_term
    procedure :: uze_theq_term
    procedure :: pak_df_obs => uze_df_obs
    procedure :: pak_rp_obs => uze_rp_obs
    procedure :: pak_bd_obs => uze_bd_obs
    procedure :: pak_set_stressperiod => uze_set_stressperiod
    procedure :: apt_ad_chk => uze_ad_chk
    procedure :: bnd_ac => uze_ac
    procedure :: bnd_mc => uze_mc
    procedure :: get_mvr_depvar
    procedure, private :: area_error

  end type GweUzeType

contains

  !> @brief Create a new UZE package
  !<
  subroutine uze_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
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
    type(GweUzeType), pointer :: uzeobj
    !
    ! -- Allocate the object and assign values to object variables
    allocate (uzeobj)
    packobj => uzeobj
    !
    ! -- Create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype)
    packobj%text = text
    !
    ! -- Allocate scalars
    call uzeobj%allocate_scalars()
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

    ! -- Store pointer to flow model interface.  When the GwfGwt exchange is
    !    created, it sets fmi%bndlist so that the GWT model has access to all
    !    the flow packages
    uzeobj%fmi => fmi
    !
    ! -- Store pointer to governing equation scale factor
    uzeobj%eqnsclfac => eqnsclfac
    !
    ! -- Store pointer to shared data module for accessing cpw, rhow
    !    for the budget calculations, and for accessing the latent heat of
    !    vaporization
    uzeobj%gwecommon => gwecommon
    !
    ! -- Set labels that will be used in generalized APT class
    uzeobj%depvartype = dvt
    uzeobj%depvarunit = dvu
    uzeobj%depvarunitabbrev = dvua
  end subroutine uze_create

  !> @brief Find corresponding uze package
  !<
  subroutine find_uze_package(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use SimVariablesModule, only: errmsg
    ! -- dummy
    class(GweUzeType) :: this
    ! -- local
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
      if (associated(this%flowbudptr)) then
        found = .true.
      end if
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
            type is (UzfType)
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
      case ('INFILTRATION')
        this%idxbudinfl = ip
        this%idxbudssm(ip) = 0
      case ('REJ-INF')
        this%idxbudrinf = ip
        this%idxbudssm(ip) = 0
      case ('UZET')
        this%idxbuduzet = ip
        this%idxbudssm(ip) = 0
      case ('REJ-INF-TO-MVR')
        this%idxbudritm = ip
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
      !
      write (this%iout, '(a, i0, " = ", a,/, a, i0)') &
        '  TERM ', ip, trim(adjustl(this%flowbudptr%budterm(ip)%flowtype)), &
        '   MAX NO. OF ENTRIES = ', this%flowbudptr%budterm(ip)%maxlist
    end do
    write (this%iout, '(a, //)') 'DONE PROCESSING '//ftype//' INFORMATION'
    !
    ! -- Thermal equilibration term
    this%idxbudtheq = this%flowbudptr%nbudterm + 1
  end subroutine find_uze_package

  !> @brief Add package connection to matrix.
  !!
  !! Overrides apt_ac to fold the UZE heat balance terms into the row
  !! corresponding to the host cell and enforce thermal equilibrium between
  !! UZE and the GWE cell.
  !<
  subroutine uze_ac(this, moffset, sparse)
    use MemoryManagerModule, only: mem_setptr
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(GweUzeType), intent(inout) :: this
    integer(I4B), intent(in) :: moffset !< current models starting position in global matrix
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    integer(I4B) :: i, ii
    integer(I4B) :: n !< index of a uze object within the complete list of uze objects
    integer(I4B) :: jj !<
    integer(I4B) :: nglo !< index of uze object in global matrix
    integer(I4B) :: jglo !< host cell's position in global matrix for a uze object
    integer(I4B) :: idxn !< used for cross-check
    integer(I4B) :: idxjj !< used for cross-check
    integer(I4B) :: idxnglo !< used for cross-check
    integer(I4B) :: idxjglo !< used for cross-check
    !
    ! -- Add package rows to sparse
    if (this%imatrows /= 0) then
      !
      ! -- Diagonal on the row assoc. with the uze feature
      do n = 1, this%ncv
        nglo = moffset + this%dis%nodes + this%ioffset + n
        call sparse%addconnection(nglo, nglo, 1)
      end do
      !
      ! -- Add uze-to-gwe connections. For uze, this particular do loop
      !    is the same as its counterpart in apt_ac.
      !    nlist: number of gwe cells with a connection to at least one uze object
      do i = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(i) !< uze object position within uze object list
        jj = this%flowbudptr%budterm(this%idxbudgwf)%id2(i) !< position of gwe cell to which uze feature is connected
        nglo = moffset + this%dis%nodes + this%ioffset + n !< uze feature position
        jglo = moffset + jj !< gwe cell position
        call sparse%addconnection(nglo, jglo, 1)
        call sparse%addconnection(jglo, nglo, 1)
      end do
      !
      ! -- For uze, add feature-to-feature connections (i.e.,
      !    vertically stacked UZ objects) to row corresponding
      !    to the host cell. Terms added to the row associated with host
      !    cell are added to columns associated with other uze features.
      !    This approach deviates from the approach taken in apt_ac.
      if (this%idxbudfjf /= 0) then
        do i = 1, this%flowbudptr%budterm(this%idxbudfjf)%maxlist
          n = this%flowbudptr%budterm(this%idxbudfjf)%id1(i) !< position of currently considered uze feature
          jj = this%flowbudptr%budterm(this%idxbudfjf)%id2(i) !< position of connected uze feature
          nglo = moffset + this%dis%nodes + this%ioffset + n !< global position of currently considered uze feature
          jglo = moffset + this%dis%nodes + this%ioffset + jj !< global position of connected uze feature
          ! -- If connected uze feature is upstream, find cell that hosts currently
          !    considered uze feature and add connection to that cell's row
          do ii = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist !< uze object id among uze objects
            idxn = this%flowbudptr%budterm(this%idxbudgwf)%id1(ii) !< uze object position within uze object list
            idxjj = this%flowbudptr%budterm(this%idxbudgwf)%id2(ii) !< position of gwe cell to which uze feature is connected
            idxnglo = moffset + this%dis%nodes + this%ioffset + idxn !< uze feature global position
            idxjglo = moffset + idxjj !< gwe cell global position
            if (nglo == idxnglo) exit
          end do
          call sparse%addconnection(idxjglo, jglo, 1)
        end do
      end if
    end if
  end subroutine uze_ac

  !> @brief Map package connection to matrix
  !<
  subroutine uze_mc(this, moffset, matrix_sln)
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(GweUzeType), intent(inout) :: this
    integer(I4B), intent(in) :: moffset
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: n, j, iglo, jglo
    integer(I4B) :: idxn, idxj, idxiglo, idxjglo
    integer(I4B) :: ipos, idxpos
    !
    ! -- Allocate memory for index arrays
    call this%apt_allocate_index_arrays()
    !
    ! -- Store index positions
    if (this%imatrows /= 0) then
      !
      ! -- Find the position of each connection in the global ia, ja structure
      !    and store them in idxglo.  idxglo allows this model to insert or
      !    retrieve values into or from the global A matrix
      !    apt rows
      !
      ! -- Feature diagonal in global matrix
      do n = 1, this%ncv
        iglo = moffset + this%dis%nodes + this%ioffset + n
        this%idxpakdiag(n) = matrix_sln%get_position_diag(iglo)
      end do
      !
      ! -- Cell to feature connection in global matrix
      do ipos = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(ipos) !< feature number
        j = this%flowbudptr%budterm(this%idxbudgwf)%id2(ipos) !< cell number
        iglo = moffset + this%dis%nodes + this%ioffset + n !< feature row index
        jglo = j + moffset !< cell row index
        ! -- Note that this is where idxlocnode is set for uze; it is set
        !    to the host cell local row index rather than the feature local
        !    row index
        this%idxlocnode(n) = j ! kluge note: do we want to introduce a new array instead of co-opting idxlocnode???
        ! -- For connection ipos in list of feature-cell connections,
        !    global positions of feature-row diagonal and off-diagonal
        !    corresponding to the cell
        this%idxdglo(ipos) = matrix_sln%get_position_diag(iglo)
        this%idxoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
      end do
      !
      ! -- Feature to cell connection in global matrix
      do ipos = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(ipos) !< feature number
        j = this%flowbudptr%budterm(this%idxbudgwf)%id2(ipos) !< cell number
        iglo = j + moffset !< cell row index
        jglo = moffset + this%dis%nodes + this%ioffset + n !< feature row index
        ! -- For connection ipos in list of feature-cell connections,
        !    global positions of cell-row diagonal and off-diagonal
        !    corresponding to the feature
        this%idxsymdglo(ipos) = matrix_sln%get_position_diag(iglo)
        this%idxsymoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
      end do
      !
      ! -- Feature to feature connection in global matrix
      if (this%idxbudfjf /= 0) then
        do ipos = 1, this%flowbudptr%budterm(this%idxbudfjf)%nlist
          n = this%flowbudptr%budterm(this%idxbudfjf)%id1(ipos) !< number of currently considered uze feature
          j = this%flowbudptr%budterm(this%idxbudfjf)%id2(ipos) !< number of connected uze feature
          iglo = moffset + this%dis%nodes + this%ioffset + n !< global position of currently considered uze feature
          jglo = moffset + this%dis%nodes + this%ioffset + j !< global position of connected uze feature
          ! -- If connected uze feature is upstream, find cell that hosts currently
          !    considered uze feature and map connection to that cell's row
          do idxpos = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
            idxn = this%flowbudptr%budterm(this%idxbudgwf)%id1(idxpos) !< feature number
            idxj = this%flowbudptr%budterm(this%idxbudgwf)%id2(idxpos) !< cell number)
            idxjglo = moffset + this%dis%nodes + this%ioffset + idxn !< feature row index
            idxiglo = moffset + idxj !< cell row index
            if (idxjglo == iglo) exit
          end do
          ! -- For connection ipos in list of feature-feature connections,
          !    global positions of host-cell-row entries corresponding to
          !    (in the same columns as) the feature-id1-row diagonal and the
          !    feature-id1-row off-diagonal corresponding to feature id2
          this%idxfjfdglo(ipos) = matrix_sln%get_position_diag(idxiglo)
          this%idxfjfoffdglo(ipos) = matrix_sln%get_position(idxiglo, jglo)
        end do
      end if
    end if
  end subroutine uze_mc

  !> @brief Add matrix terms related to UZE
  !!
  !! This will be called from TspAptType%apt_fc_expanded()
  !! in order to add matrix terms specifically for this package
  !<
  subroutine uze_fc_expanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(GweUzeType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, n, n1, n2
    integer(I4B) :: iloc
    integer(I4B) :: iposd, iposoffd
    integer(I4B) :: ipossymoffd
    real(DP) :: cold
    real(DP) :: qbnd
    real(DP) :: omega
    real(DP) :: rrate
    real(DP) :: rhsval
    real(DP) :: hcofval
    !
    ! -- Add infiltration contribution
    !    uze does not put feature balance coefficients in the row
    !    associated with the feature.  Instead, these coefficients are
    !    moved into the row associated with cell hosting the uze feature
    if (this%idxbudinfl /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudinfl)%nlist
        call this%uze_infl_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1) !< for uze idxlocnode stores the host cell local row index
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymoffd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add rejected infiltration contribution
    if (this%idxbudrinf /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrinf)%nlist
        call this%uze_rinf_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1) ! for uze idxlocnode stores the host cell local row index
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymoffd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add unsaturated et contribution
    if (this%idxbuduzet /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbuduzet)%nlist
        call this%uze_uzet_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1) ! for uze idxlocnode stores the host cell local row index
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymoffd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add rejected infiltration to mover contribution
    if (this%idxbudritm /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudritm)%nlist
        call this%uze_ritm_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1) ! for uze idxlocnode stores the host cell local row index
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymoffd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- For UZE, content of apt_fc_expanded placed here as the approach is to
    !    completely override apt_fc_expanded() with what follows
    !
    ! -- Mass (or energy) storage in features
    do n = 1, this%ncv
      cold = this%xoldpak(n)
      iloc = this%idxlocnode(n) ! for uze idxlocnode stores the host cell local row index
      ipossymoffd = this%idxsymoffdglo(n)
      call this%apt_stor_term(n, n1, n2, rrate, rhsval, hcofval)
      call matrix_sln%add_value_pos(ipossymoffd, hcofval)
      rhs(iloc) = rhs(iloc) + rhsval
    end do
    !
    ! -- Add to mover contribution
    if (this%idxbudtmvr /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudtmvr)%nlist
        call this%apt_tmvr_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1) ! for uze, idxlocnode stores the host cell local row index
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymoffd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- Add from mover contribution
    if (this%idxbudfmvr /= 0) then
      do n = 1, this%ncv
        rhsval = this%qmfrommvr(n) ! kluge note: presumably already in terms of energy
        iloc = this%idxlocnode(n) ! for uze idxlocnode stores the host cell local row index
        rhs(iloc) = rhs(iloc) - rhsval
      end do
    end if
    !
    ! -- Go through each apt-gwf connection
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      !
      ! -- Set n to feature number and process if active feature
      n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      if (this%iboundpak(n) /= 0) then
        !
        ! -- This code altered from its counterpart appearing in apt; this equates
        !    uze temperature to cell temperature using the feature's row
        iposd = this%idxdglo(j)
        iposoffd = this%idxoffdglo(j)
        call matrix_sln%add_value_pos(iposd, DONE)
        call matrix_sln%add_value_pos(iposoffd, -DONE)
      end if
    end do
    !
    ! -- Go through each apt-apt connection
    if (this%idxbudfjf /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfjf)%nlist
        n1 = this%flowbudptr%budterm(this%idxbudfjf)%id1(j)
        n2 = this%flowbudptr%budterm(this%idxbudfjf)%id2(j)
        qbnd = this%flowbudptr%budterm(this%idxbudfjf)%flow(j)
        if (qbnd <= DZERO) then
          omega = DONE
        else
          omega = DZERO
        end if
        iposd = this%idxfjfdglo(j) !< position of feature-id1 column in feature id1's host-cell row
        iposoffd = this%idxfjfoffdglo(j) !< position of feature-id2 column in feature id1's host-cell row
        call matrix_sln%add_value_pos(iposd, omega * qbnd * this%eqnsclfac)
        call matrix_sln%add_value_pos(iposoffd, &
                                      (DONE - omega) * qbnd * this%eqnsclfac)
      end do
    end if
  end subroutine uze_fc_expanded

  !> @brief Explicit solve
  !!
  !! There should be no explicit solve for uze.  However, if there were, then
  !! this subroutine would add terms specific to the unsaturated zone to the
  !! explicit unsaturated-zone solve
  subroutine uze_solve(this)
    ! -- dummy
    class(GweUzeType) :: this
    ! -- local
    integer(I4B) :: j
    integer(I4B) :: n1, n2
    real(DP) :: rrate
    !
    ! -- Add infiltration contribution
    if (this%idxbudinfl /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudinfl)%nlist
        call this%uze_infl_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add rejected infiltration contribution
    if (this%idxbudrinf /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudrinf)%nlist
        call this%uze_rinf_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add unsaturated et contribution
    if (this%idxbuduzet /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbuduzet)%nlist
        call this%uze_uzet_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- Add rejected infiltration to mover contribution
    if (this%idxbudritm /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudritm)%nlist
        call this%uze_ritm_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
  end subroutine uze_solve

  !> @brief Return the number of UZE-specific budget terms
  !!
  !! Function to return the number of budget terms just for this package.
  !! This overrides function in parent.
  !<
  function uze_get_nbudterms(this) result(nbudterms)
    ! -- dummy
    class(GweUzeType) :: this
    ! -- Return
    integer(I4B) :: nbudterms
    !
    ! -- Number of budget terms is 5
    nbudterms = 0
    if (this%idxbudinfl /= 0) nbudterms = nbudterms + 1
    if (this%idxbudrinf /= 0) nbudterms = nbudterms + 1
    if (this%idxbuduzet /= 0) nbudterms = nbudterms + 1
    if (this%idxbudritm /= 0) nbudterms = nbudterms + 1
    if (this%idxbudtheq /= 0) nbudterms = nbudterms + 1
  end function uze_get_nbudterms

  !> @brief Override similarly named function in APT
  !!
  !! Set the temperature to be used by MVE as the user-specified
  !! temperature applied to the infiltration
  !<
  function get_mvr_depvar(this)
    ! -- dummy
    class(GweUzeType) :: this
    ! -- return
    real(dp), dimension(:), contiguous, pointer :: get_mvr_depvar
    !
    get_mvr_depvar => this%tempinfl
  end function get_mvr_depvar

  !> @brief Setup budget object
  !!
  !! Set up the budget object that stores all the unsaturated-zone
  !! flows
  !<
  subroutine uze_setup_budobj(this, idx)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(GweUzeType) :: this
    integer(I4B), intent(inout) :: idx
    ! -- local
    integer(I4B) :: maxlist, naux
    character(len=LENBUDTXT) :: text
    !
    ! -- Infiltration
    text = '    INFILTRATION'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudinfl)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
    !
    ! -- Rejected infiltration (Hortonian flow)
    if (this%idxbudrinf /= 0) then
      text = '         REJ-INF'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudrinf)%maxlist
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
    ! -- Evapotranspiration from the unsaturated zone
    if (this%idxbuduzet /= 0) then
      text = '            UZET'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbuduzet)%maxlist
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
    ! -- Rejected infiltration that is subsequently transferred by MVR
    if (this%idxbudritm /= 0) then
      text = '  INF-REJ-TO-MVR'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudritm)%maxlist
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
    ! -- Energy transferred to solid phase by the thermal equilibrium assumption
    text = '   THERMAL-EQUIL'
    idx = idx + 1
    ! -- use dimension of GWF term
    maxlist = this%flowbudptr%budterm(this%idxbudgwf)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux)
  end subroutine uze_setup_budobj

  !> @brief Fill UZE budget object
  !!
  !! Copy flow terms into this%budobj
  !<
  subroutine uze_fill_budobj(this, idx, x, flowja, ccratin, ccratout)
    ! -- modules
    ! -- dummy
    class(GweUzeType) :: this
    integer(I4B), intent(inout) :: idx
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    real(DP), intent(inout) :: ccratin
    real(DP), intent(inout) :: ccratout
    ! -- local
    integer(I4B) :: j, n1, n2, indx
    integer(I4B) :: nlist, nlen
    integer(I4B) :: igwfnode
    integer(I4B) :: idiag
    real(DP) :: q
    real(DP), dimension(:), allocatable :: budresid
    !
    allocate (budresid(this%ncv))
    do n1 = 1, this%ncv
      budresid(n1) = DZERO
    end do
    !
    indx = 0
    !
    ! -- FLOW JA FACE into budresid
    nlen = 0
    if (this%idxbudfjf /= 0) then
      nlen = this%flowbudptr%budterm(this%idxbudfjf)%maxlist
    end if
    if (nlen > 0) then
      indx = indx + 1
      nlist = this%budobj%budterm(indx)%nlist
      do j = 1, nlist
        n1 = this%budobj%budterm(indx)%id1(j)
        n2 = this%budobj%budterm(indx)%id2(j)
        if (n1 < n2) then
          q = this%budobj%budterm(indx)%flow(j)
          budresid(n1) = budresid(n1) + q
          budresid(n2) = budresid(n2) - q
        end if
      end do
    end if
    !
    ! -- GWF (LEAKAGE) into budresid
    indx = indx + 1
    nlist = this%budobj%budterm(indx)%nlist
    do j = 1, nlist
      n1 = this%budobj%budterm(indx)%id1(j)
      q = this%budobj%budterm(indx)%flow(j)
      budresid(n1) = budresid(n1) + q
    end do
    !
    ! -- Skip individual package terms
    indx = this%idxlastpak
    !
    ! -- STORAGE into budresid
    indx = indx + 1
    do n1 = 1, this%ncv
      q = this%budobj%budterm(indx)%flow(n1)
      budresid(n1) = budresid(n1) + q
    end do
    !
    ! -- TO MOVER into budresid
    if (this%idxbudtmvr /= 0) then
      indx = indx + 1
      nlist = this%budobj%budterm(indx)%nlist
      do j = 1, nlist
        n1 = this%budobj%budterm(indx)%id1(j)
        q = this%budobj%budterm(indx)%flow(j)
        budresid(n1) = budresid(n1) + q
      end do
    end if
    !
    ! -- FROM MOVER into budresid
    if (this%idxbudfmvr /= 0) then
      indx = indx + 1
      nlist = this%budobj%budterm(indx)%nlist
      do j = 1, nlist
        n1 = this%budobj%budterm(indx)%id1(j)
        q = this%budobj%budterm(indx)%flow(j)
        budresid(n1) = budresid(n1) + q
      end do
    end if
    !
    ! -- CONSTANT FLOW into budresid
    indx = indx + 1
    do n1 = 1, this%ncv
      q = this%budobj%budterm(indx)%flow(n1)
      budresid(n1) = budresid(n1) + q
    end do
    !
    ! -- AUXILIARY VARIABLES into budresid
    ! -- (No flows associated with these)
    !
    ! -- Individual package terms processed last
    !
    ! -- Infiltration
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudinfl)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      call this%uze_infl_term(j, n1, n2, q)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      budresid(n1) = budresid(n1) + q
    end do
    !
    ! -- Rej-Inf
    if (this%idxbudrinf /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudrinf)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%uze_rinf_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
        budresid(n1) = budresid(n1) + q
      end do
    end if
    !
    ! -- UZET
    if (this%idxbuduzet /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbuduzet)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%uze_uzet_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
        budresid(n1) = budresid(n1) + q
      end do
    end if
    !
    ! -- Rej-Inf-To-MVR
    if (this%idxbudritm /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudritm)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%uze_ritm_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
        budresid(n1) = budresid(n1) + q
      end do
    end if
    !
    ! -- Thermal-Equil
    ! -- processed last because it is calculated from the residual
    idx = idx + 1
    nlist = this%flowbudptr%budterm(this%idxbudgwf)%nlist
    call this%budobj%budterm(idx)%reset(nlist)
    do j = 1, nlist
      n1 = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(j)
      q = -budresid(n1)
      call this%uze_theq_term(j, n1, igwfnode, q)
      call this%budobj%budterm(idx)%update_term(n1, igwfnode, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      if (this%iboundpak(n1) /= 0) then
        ! -- Contribution to gwe cell budget
        this%simvals(n1) = this%simvals(n1) - q
        idiag = this%dis%con%ia(igwfnode)
        flowja(idiag) = flowja(idiag) - q
      end if
    end do
    !
    deallocate (budresid)
  end subroutine uze_fill_budobj

  !> @brief Allocate scalars
  !!
  !! Allocate scalars specific to UZE package
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweUzeType) :: this
    ! -- local
    !
    ! -- Allocate scalars in TspAptType
    call this%TspAptType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%idxbudinfl, 'IDXBUDINFL', this%memoryPath)
    call mem_allocate(this%idxbudrinf, 'IDXBUDRINF', this%memoryPath)
    call mem_allocate(this%idxbuduzet, 'IDXBUDUZET', this%memoryPath)
    call mem_allocate(this%idxbudritm, 'IDXBUDRITM', this%memoryPath)
    call mem_allocate(this%idxbudtheq, 'IDXBUDTHEQ', this%memoryPath)
    !
    ! -- Initialize
    this%idxbudinfl = 0
    this%idxbudrinf = 0
    this%idxbuduzet = 0
    this%idxbudritm = 0
    this%idxbudtheq = 0
  end subroutine allocate_scalars

  !> @brief Allocate arrays
  !!
  !! Allocate arrays used by the UZE package
  !<
  subroutine uze_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GweUzeType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- Time series
    call mem_allocate(this%tempinfl, this%ncv, 'TEMPINFL', this%memoryPath)
    call mem_allocate(this%tempuzet, this%ncv, 'TEMPUZET', this%memoryPath)
    !
    ! -- Call standard TspAptType allocate arrays
    call this%TspAptType%apt_allocate_arrays()
    !
    ! -- Initialize
    do n = 1, this%ncv
      this%tempinfl(n) = DZERO
      this%tempuzet(n) = DZERO
    end do
  end subroutine uze_allocate_arrays

  !> @brief Deallocate memory
  !<
  subroutine uze_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GweUzeType) :: this
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%idxbudinfl)
    call mem_deallocate(this%idxbudrinf)
    call mem_deallocate(this%idxbuduzet)
    call mem_deallocate(this%idxbudritm)
    call mem_deallocate(this%idxbudtheq)
    !
    ! -- Deallocate time series
    call mem_deallocate(this%tempinfl)
    call mem_deallocate(this%tempuzet)
    !
    ! -- Deallocate scalars in TspAptType
    call this%TspAptType%bnd_da()
  end subroutine uze_da

  !> @brief Infiltration term
  !!
  !! Accounts for energy added to the subsurface via infiltration, for example,
  !! energy entering the model domain via rainfall or irrigation.
  !<
  subroutine uze_infl_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweUzeType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudinfl)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudinfl)%id2(ientry)
    !
    ! -- Note that qbnd is negative for negative infiltration
    qbnd = this%flowbudptr%budterm(this%idxbudinfl)%flow(ientry)
    if (qbnd < DZERO) then
      ctmp = this%xnewpak(n1)
      h = qbnd
      r = DZERO
    else
      ctmp = this%tempinfl(n1)
      h = DZERO
      r = -qbnd * ctmp
    end if
    if (present(rrate)) rrate = qbnd * ctmp * this%eqnsclfac
    if (present(rhsval)) rhsval = r * this%eqnsclfac
    if (present(hcofval)) hcofval = h * this%eqnsclfac
  end subroutine uze_infl_term

  !> @brief Rejected infiltration term
  !!
  !! Accounts for energy that is added to the model from specifying an
  !! infiltration rate and temperature, but is subsequently removed from
  !! the model as that portion of the infiltration that is rejected (and
  !! NOT transferred to another advanced package via the MVR/MVT packages).
  !<
  subroutine uze_rinf_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweUzeType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudrinf)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudrinf)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudrinf)%flow(ientry)
    ctmp = this%tempinfl(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine uze_rinf_term

  !> @brief Evapotranspiration from the unsaturated-zone term
  !!
  !! Accounts for thermal cooling in the unsaturated zone as a result of
  !! evapotranspiration from the unsaturated zone.  Amount of water converted
  !! to vapor phase (UZET) determined by GWF model
  !<
  subroutine uze_uzet_term(this, ientry, n1, n2, rrate, rhsval, hcofval)
    ! -- dummy
    class(GweUzeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    ! -- local
    real(DP) :: qbnd
    real(DP) :: ctmp
    real(DP) :: omega
    !
    n1 = this%flowbudptr%budterm(this%idxbuduzet)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbuduzet)%id2(ientry)
    ! -- Note that qbnd is negative for uzet
    qbnd = this%flowbudptr%budterm(this%idxbuduzet)%flow(ientry)
    ctmp = this%tempuzet(n1)
    if (this%xnewpak(n1) < ctmp) then
      omega = DONE
    else
      omega = DZERO
    end if
    if (present(rrate)) &
      rrate = (omega * qbnd * this%xnewpak(n1) + &
               (DONE - omega) * qbnd * ctmp) * this%eqnsclfac
    if (present(rhsval)) rhsval = -(DONE - omega) * qbnd * ctmp * this%eqnsclfac
    if (present(hcofval)) hcofval = omega * qbnd * this%eqnsclfac
  end subroutine uze_uzet_term

  !> @brief Rejected infiltration to MVR/MVT term
  !!
  !! Accounts for energy that is added to the model from specifying an
  !! infiltration rate and temperature, but does not infiltrate into the
  !! subsurface.  This subroutine is called when the rejected infiltration
  !! is transferred to another advanced package via the MVR/MVT packages.
  !<
  subroutine uze_ritm_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- dummy
    class(GweUzeType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudritm)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudritm)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudritm)%flow(ientry)
    ctmp = this%tempinfl(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine uze_ritm_term

  !> @brief Heat transferred through thermal equilibrium with the solid phase
  !!
  !! Accounts for the transfer of energy from the liquid phase to the solid
  !! phase as a result of the instantaneous thermal equilibrium assumption.
  !<
  subroutine uze_theq_term(this, ientry, n1, n2, rrate)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(GweUzeType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout) :: rrate
    ! -- local
    real(DP) :: r
    integer(I4B) :: i
    character(len=LENBUDTXT) :: flowtype
    !
    r = DZERO
    n1 = this%flowbudptr%budterm(this%idxbudgwf)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudgwf)%id2(ientry)
    if (this%iboundpak(n1) /= 0) then
      do i = 1, this%budobj%nbudterm
        flowtype = this%budobj%budterm(i)%flowtype
        select case (trim(adjustl(flowtype)))
        case ('THERMAL-EQUIL')
          ! -- Skip
          continue
        case default
          r = r - this%budobj%budterm(i)%flow(ientry)
        end select
      end do
    end if
    rrate = r
  end subroutine uze_theq_term

  !> @brief Define UZE Observation
  !!
  !! This subroutine:
  !!   - Stores observation types supported by the parent APT package.
  !!   - Overrides BndType%bnd_df_obs
  !<
  subroutine uze_df_obs(this)
    ! -- dummy
    class(GweUzeType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Store obs type and assign procedure pointer
    !    for temperature observation type.
    call this%obs%StoreObsType('temperature', .false., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for flow between uze cells.
    call this%obs%StoreObsType('flow-ja-face', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID12
    !
    ! -- Store obs type and assign procedure pointer
    !    for from-mvr observation type.
    call this%obs%StoreObsType('from-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- to-mvr not supported for uze
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
    !    for observation type: uze
    call this%obs%StoreObsType('uze', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('infiltration', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('rej-inf', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('uzet', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('rej-inf-to-mvr', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
    !
    ! -- Store obs type and assign procedure pointer
    !    for observation type.
    call this%obs%StoreObsType('thermal-equil', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => apt_process_obsID
  end subroutine uze_df_obs

  !> @brief Process package specific obs
  !!
  !! Method to process specific observations for this package.
  !<
  subroutine uze_rp_obs(this, obsrv, found)
    ! -- dummy
    class(GweUzeType), intent(inout) :: this !< package class
    type(ObserveType), intent(inout) :: obsrv !< observation object
    logical, intent(inout) :: found !< indicate whether observation was found
    !
    found = .true.
    select case (obsrv%ObsTypeId)
    case ('INFILTRATION')
      call this%rp_obs_byfeature(obsrv)
    case ('REJ-INF')
      call this%rp_obs_byfeature(obsrv)
    case ('UZET')
      call this%rp_obs_byfeature(obsrv)
    case ('REJ-INF-TO-MVR')
      call this%rp_obs_byfeature(obsrv)
    case ('THERMAL-EQUIL')
      call this%rp_obs_byfeature(obsrv)
    case default
      found = .false.
    end select
  end subroutine uze_rp_obs

  !> @brief Calculate observation value and pass it back to APT
  !<
  subroutine uze_bd_obs(this, obstypeid, jj, v, found)
    ! -- dummy
    class(GweUzeType), intent(inout) :: this
    character(len=*), intent(in) :: obstypeid
    real(DP), intent(inout) :: v
    integer(I4B), intent(in) :: jj
    logical, intent(inout) :: found
    ! -- local
    integer(I4B) :: n1, n2
    !
    found = .true.
    select case (obstypeid)
    case ('INFILTRATION')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudinfl > 0) then
        call this%uze_infl_term(jj, n1, n2, v)
      end if
    case ('REJ-INF')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudrinf > 0) then
        call this%uze_rinf_term(jj, n1, n2, v)
      end if
    case ('UZET')
      if (this%iboundpak(jj) /= 0 .and. this%idxbuduzet > 0) then
        call this%uze_uzet_term(jj, n1, n2, v)
      end if
    case ('REJ-INF-TO-MVR')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudritm > 0) then
        call this%uze_ritm_term(jj, n1, n2, v)
      end if
    case ('THERMAL-EQUIL')
      if (this%iboundpak(jj) /= 0 .and. this%idxbudtheq > 0) then
        call this%uze_theq_term(jj, n1, n2, v)
      end if
    case default
      found = .false.
    end select
  end subroutine uze_bd_obs

  !> @brief Check if UZF object area is not equal to the cell area
  !!
  !! GWE equilibrates the temperature of a UZE object with the host cell. UZE
  !! assumes that the area associated with a specific UZE object is equal to
  !! the area of the host cell. When this condition is not true, the code
  !! exits with an appropriate message.
  !<
  subroutine uze_ad_chk(this)
    ! -- modules
    use ConstantsModule, only: IZERO
    use MathUtilModule, only: is_close
    use SimModule, only: count_errors
    ! -- dummy
    class(GweUzeType), intent(inout) :: this
    ! -- local
    integer(I4B) :: nuz
    integer(I4B) :: n
    integer(I4B) :: igwfnode
    real(DP) :: carea
    real(DP) :: uzarea

    nuz = this%flowbudptr%budterm(this%idxbudgwf)%maxlist

    ! cycle through uze objects, stop at first occurrence of more than one
    ! uze object in a cell
    do n = 1, nuz
      igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(n)
      carea = this%dis%area(igwfnode)
      uzarea = this%flowbudptr%budterm(this%idxbudgwf)%auxvar(1, n)
      ! compare areas
      if (.not. is_close(carea, uzarea)) then
        call this%area_error(igwfnode)
      end if
    end do
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine uze_ad_chk

  !> @brief Print and store error msg indicating area of UZF object is not
  !! equal to that of the host cell
  !<
  subroutine area_error(this, iloc)
    ! -- modules
    use SimVariablesModule, only: errmsg
    ! -- dummy
    class(GweUzeType) :: this
    integer(I4B) :: iloc
    ! -- local
    character(len=30) :: nodestr
    !
    call this%dis%noder_to_string(iloc, nodestr)
    write (errmsg, '(3a)') &
      'In a GWE model, the area of every UZF object must be equal to that of &
      &the host cell. This condition is violated in cell ', &
       trim(adjustl(nodestr)), '. Check use of AUXMULTNAME option in UZF &
      &package.'
    call store_error(errmsg)
  end subroutine area_error

  !> @brief Sets the stress period attributes for keyword use.
  !<
  subroutine uze_set_stressperiod(this, itemno, keyword, found)
    ! -- modules
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(GweUzeType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    character(len=*), intent(in) :: keyword
    logical, intent(inout) :: found
    ! -- local
    character(len=LINELENGTH) :: temp_text
    integer(I4B) :: ierr
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    !
    ! INFILTRATION <infiltration>
    ! UZET <uzet>
    !
    found = .true.
    select case (keyword)
    case ('INFILTRATION')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(temp_text)
      jj = 1
      bndElem => this%tempinfl(itemno)
      call read_value_or_time_series_adv(temp_text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'INFILTRATION')
    case ('UZET')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(temp_text)
      jj = 1
      bndElem => this%tempuzet(itemno)
      call read_value_or_time_series_adv(temp_text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, 'UZET')
    case default
      !
      ! -- Keyword not recognized so return to caller with found = .false.
      found = .false.
    end select
    !
999 continue
  end subroutine uze_set_stressperiod

end module GweUzeModule
