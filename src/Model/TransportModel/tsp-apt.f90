! -- Advanced Package Transport Module
! -- This module contains most of the routines for simulating transport
! -- through the advanced packages.
! -- Future work:
!      * support decay, sorption
!      * dispersion in SFT and UZT?
!
! AFP flows (flowbudptr)    index var     ATP term              Transport Type
!---------------------------------------------------------------------------------

! -- specialized terms in the flow budget
! FLOW-JA-FACE              idxbudfjf     FLOW-JA-FACE          cv2cv
! GWF (aux FLOW-AREA)       idxbudgwf     GWF                   cv2gwf
! STORAGE (aux VOLUME)      idxbudsto     none                  used for cv volumes
! FROM-MVR                  idxbudfmvr    FROM-MVR              q * cext = this%qfrommvr(:)   ! rhow*cpw is applied to various terms for heat transport
! TO-MVR                    idxbudtmvr    TO-MVR                q * cfeat

! -- generalized source/sink terms (except ET?)
! RAINFALL                  idxbudrain    RAINFALL              q * crain
! EVAPORATION               idxbudevap    EVAPORATION           cfeat<cevap: q*cfeat, else: q*cevap  ! latent heat may be applied for evaporative cooling
! RUNOFF                    idxbudroff    RUNOFF                q * croff
! EXT-INFLOW                idxbudiflw    EXT-INFLOW            q * ciflw
! WITHDRAWAL                idxbudwdrl    WITHDRAWAL            q * cfeat
! EXT-OUTFLOW               idxbudoutf    EXT-OUTFLOW           q * cfeat

! -- terms from a flow file that should be skipped
! CONSTANT                  none          none                  none
! AUXILIARY                 none          none                  none

! -- terms that are written to the transport budget file
! none                      none          STORAGE (aux MASS)    dM/dt
! none                      none          AUXILIARY             none
! none                      none          CONSTANT              accumulate
!
!
module TspAptModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, DEP20, LENFTYPE, LINELENGTH, &
                             LENBOUNDNAME, LENPACKAGENAME, NAMEDBOUNDFLAG, &
                             DNODATA, TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL, &
                             LENAUXNAME, LENVARNAME
  use SimModule, only: store_error, store_error_unit, count_errors
  use SimVariablesModule, only: errmsg
  use BndModule, only: BndType
  use TspFmiModule, only: TspFmiType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use BudgetTermModule, only: BudgetTermType
  use TableModule, only: TableType, table_cr
  use ObserveModule, only: ObserveType
  use InputOutputModule, only: extract_idnum_or_bndname, str_pad_left
  use BaseDisModule, only: DisBaseType
  use MatrixBaseModule

  implicit none

  public :: TspAptType
  public :: apt_process_obsID
  public :: apt_process_obsID12

  character(len=LENFTYPE) :: ftype = 'APT'
  character(len=LENVARNAME) :: text = '             APT'

  type, extends(BndType) :: TspAptType

    character(len=LENPACKAGENAME) :: flowpackagename = '' !< name of corresponding flow package
    character(len=8), &
      dimension(:), pointer, contiguous :: status => null() !< active, inactive, constant
    character(len=LENAUXNAME) :: cauxfpconc = '' !< name of aux column in flow package auxvar array for concentration (or temperature)
    integer(I4B), pointer :: iauxfpconc => null() !< column in flow package bound array to insert concs
    integer(I4B), pointer :: imatrows => null() !< if active, add new rows to matrix
    integer(I4B), pointer :: iprconc => null() !< print conc to listing file
    integer(I4B), pointer :: iconcout => null() !< unit number for conc output file
    integer(I4B), pointer :: ibudgetout => null() !< unit number for budget output file
    integer(I4B), pointer :: ibudcsv => null() !< unit number for csv budget output file
    integer(I4B), pointer :: ncv => null() !< number of control volumes
    integer(I4B), pointer :: igwfaptpak => null() !< package number of corresponding this package
    integer(I4B), pointer :: idxprepak => null() !< budget-object index that precedes package-specific budget objects
    integer(I4B), pointer :: idxlastpak => null() !< budget-object index of last package-specific budget object
    real(DP), dimension(:), pointer, contiguous :: strt => null() !< starting feature concentration (or temperature)
    integer(I4B), dimension(:), pointer, contiguous :: idxlocnode => null() !< map position in global rhs and x array of pack entry
    integer(I4B), dimension(:), pointer, contiguous :: idxpakdiag => null() !< map diag position of feature in global amat
    integer(I4B), dimension(:), pointer, contiguous :: idxdglo => null() !< map position in global array of package diagonal row entries
    integer(I4B), dimension(:), pointer, contiguous :: idxoffdglo => null() !< map position in global array of package off diagonal row entries
    integer(I4B), dimension(:), pointer, contiguous :: idxsymdglo => null() !< map position in global array of package diagonal entries to model rows
    integer(I4B), dimension(:), pointer, contiguous :: idxsymoffdglo => null() !< map position in global array of package off diagonal entries to model rows
    integer(I4B), dimension(:), pointer, contiguous :: idxfjfdglo => null() !< map diagonal feature to feature in global amat
    integer(I4B), dimension(:), pointer, contiguous :: idxfjfoffdglo => null() !< map off diagonal feature to feature in global amat
    integer(I4B), dimension(:), pointer, contiguous :: iboundpak => null() !< package ibound
    real(DP), dimension(:), pointer, contiguous :: xnewpak => null() !< feature concentration (or temperature) for current time step
    real(DP), dimension(:), pointer, contiguous :: xoldpak => null() !< feature concentration (or temperature) from previous time step
    real(DP), dimension(:), pointer, contiguous :: dbuff => null() !< temporary storage array
    character(len=LENBOUNDNAME), &
      dimension(:), pointer, contiguous :: featname => null()
    real(DP), dimension(:), pointer, contiguous :: concfeat => null() !< concentration (or temperature) of the feature
    real(DP), dimension(:, :), pointer, contiguous :: lauxvar => null() !< auxiliary variable
    type(TspFmiType), pointer :: fmi => null() !< pointer to fmi object
    real(DP), dimension(:), pointer, contiguous :: qsto => null() !< mass (or energy) flux due to storage change
    real(DP), dimension(:), pointer, contiguous :: ccterm => null() !< mass (or energy) flux required to maintain constant concentration (or temperature)
    integer(I4B), pointer :: idxbudfjf => null() !< index of flow ja face in flowbudptr
    integer(I4B), pointer :: idxbudgwf => null() !< index of gwf terms in flowbudptr
    integer(I4B), pointer :: idxbudsto => null() !< index of storage terms in flowbudptr
    integer(I4B), pointer :: idxbudtmvr => null() !< index of to mover terms in flowbudptr
    integer(I4B), pointer :: idxbudfmvr => null() !< index of from mover terms in flowbudptr
    integer(I4B), pointer :: idxbudaux => null() !< index of auxiliary terms in flowbudptr
    integer(I4B), dimension(:), pointer, contiguous :: idxbudssm => null() !< flag that flowbudptr%buditem is a general solute source/sink
    integer(I4B), pointer :: nconcbudssm => null() !< number of concbudssm terms (columns)
    real(DP), dimension(:, :), pointer, contiguous :: concbudssm => null() !< user specified concentrations (or temperatures) for flow terms
    real(DP), dimension(:), pointer, contiguous :: qmfrommvr => null() !< a mass or energy flow coming from the mover that needs to be added
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy
    character(len=LENVARNAME) :: depvartype = '' !< stores string identifying dependent variable type, depending on model type
    character(len=LENVARNAME) :: depvarunit = '' !< "mass" or "energy"
    character(len=LENVARNAME) :: depvarunitabbrev = '' !< "M" or "E"
    !
    ! -- pointer to flow package boundary
    type(BndType), pointer :: flowpackagebnd => null()
    !
    ! -- budget objects
    type(BudgetObjectType), pointer :: budobj => null() !< apt solute budget object
    type(BudgetObjectType), pointer :: flowbudptr => null() !< GWF flow budget object
    !
    ! -- table objects
    type(TableType), pointer :: dvtab => null()

  contains

    procedure :: set_pointers => apt_set_pointers
    procedure :: bnd_ac => apt_ac
    procedure :: bnd_mc => apt_mc
    procedure :: bnd_ar => apt_ar
    procedure :: bnd_rp => apt_rp
    procedure :: bnd_ad => apt_ad
    procedure :: bnd_reset => apt_reset
    procedure :: bnd_fc => apt_fc
    procedure, public :: apt_fc_expanded ! Made public for uze
    procedure, public :: apt_ad_chk
    procedure :: pak_fc_expanded
    procedure, private :: apt_fc_nonexpanded
    procedure, public :: apt_cfupdate ! Made public for uze
    procedure :: apt_check_valid
    procedure :: apt_set_stressperiod
    procedure :: pak_set_stressperiod
    procedure :: apt_accumulate_ccterm
    procedure :: bnd_cq => apt_cq
    procedure :: bnd_ot_package_flows => apt_ot_package_flows
    procedure :: bnd_ot_dv => apt_ot_dv
    procedure :: bnd_ot_bdsummary => apt_ot_bdsummary
    procedure :: bnd_da => apt_da
    procedure :: allocate_scalars
    procedure :: apt_allocate_index_arrays
    procedure :: apt_allocate_arrays
    procedure :: find_apt_package
    procedure :: apt_solve
    procedure :: pak_solve
    procedure :: bnd_options => apt_options
    procedure :: read_dimensions => apt_read_dimensions
    procedure :: apt_read_cvs
    procedure :: read_initial_attr => apt_read_initial_attr
    procedure :: define_listlabel
    ! -- methods for observations
    procedure :: bnd_obs_supported => apt_obs_supported
    procedure :: bnd_df_obs => apt_df_obs
    procedure :: pak_df_obs
    procedure :: pak_rp_obs
    procedure :: bnd_rp_obs => apt_rp_obs
    procedure :: rp_obs_byfeature
    procedure :: rp_obs_budterm
    procedure :: rp_obs_flowjaface
    procedure :: bnd_bd_obs => apt_bd_obs
    procedure :: pak_bd_obs
    procedure :: pak_get_nbudterms
    procedure :: apt_setup_budobj
    procedure :: pak_setup_budobj
    procedure :: apt_fill_budobj
    procedure :: pak_fill_budobj
    procedure, public :: apt_get_volumes !< made public for sft/sfe
    procedure, public :: apt_stor_term
    procedure, public :: apt_tmvr_term
    procedure, public :: apt_fmvr_term !< made public for uze
    procedure, public :: apt_fjf_term !< made public for uze
    procedure, private :: apt_copy2flowp
    procedure, private :: apt_setup_tableobj
    procedure, public :: get_mvr_depvar

  end type TspAptType

contains

  !> @brief Add package connection to matrix
  !<
  subroutine apt_ac(this, moffset, sparse)
    use MemoryManagerModule, only: mem_setptr
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    integer(I4B), intent(in) :: moffset
    type(sparsematrix), intent(inout) :: sparse
    ! -- local
    integer(I4B) :: i, n
    integer(I4B) :: jj, jglo
    integer(I4B) :: nglo
    ! -- format
    !
    ! -- Add package rows to sparse
    if (this%imatrows /= 0) then
      !
      ! -- diagonal
      do n = 1, this%ncv
        nglo = moffset + this%dis%nodes + this%ioffset + n
        call sparse%addconnection(nglo, nglo, 1)
      end do
      !
      ! -- apt-gwf connections
      do i = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(i)
        jj = this%flowbudptr%budterm(this%idxbudgwf)%id2(i)
        nglo = moffset + this%dis%nodes + this%ioffset + n
        jglo = jj + moffset
        call sparse%addconnection(nglo, jglo, 1)
        call sparse%addconnection(jglo, nglo, 1)
      end do
      !
      ! -- apt-apt connections
      if (this%idxbudfjf /= 0) then
        do i = 1, this%flowbudptr%budterm(this%idxbudfjf)%maxlist
          n = this%flowbudptr%budterm(this%idxbudfjf)%id1(i)
          jj = this%flowbudptr%budterm(this%idxbudfjf)%id2(i)
          nglo = moffset + this%dis%nodes + this%ioffset + n
          jglo = moffset + this%dis%nodes + this%ioffset + jj
          call sparse%addconnection(nglo, jglo, 1)
        end do
      end if
    end if
  end subroutine apt_ac

  !> @brief Advanced package transport map package connections to matrix
  !<
  subroutine apt_mc(this, moffset, matrix_sln)
    use SparseModule, only: sparsematrix
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    integer(I4B), intent(in) :: moffset
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: n, j, iglo, jglo
    integer(I4B) :: ipos
    ! -- format
    !
    ! -- allocate memory for index arrays
    call this%apt_allocate_index_arrays()
    !
    ! -- store index positions
    if (this%imatrows /= 0) then
      !
      ! -- Find the position of each connection in the global ia, ja structure
      !    and store them in idxglo.  idxglo allows this model to insert or
      !    retrieve values into or from the global A matrix
      ! -- apt rows
      do n = 1, this%ncv
        this%idxlocnode(n) = this%dis%nodes + this%ioffset + n
        iglo = moffset + this%dis%nodes + this%ioffset + n
        this%idxpakdiag(n) = matrix_sln%get_position_diag(iglo)
      end do
      do ipos = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(ipos)
        j = this%flowbudptr%budterm(this%idxbudgwf)%id2(ipos)
        iglo = moffset + this%dis%nodes + this%ioffset + n
        jglo = j + moffset
        this%idxdglo(ipos) = matrix_sln%get_position_diag(iglo)
        this%idxoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
      end do
      !
      ! -- apt contributions to gwf portion of global matrix
      do ipos = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(ipos)
        j = this%flowbudptr%budterm(this%idxbudgwf)%id2(ipos)
        iglo = j + moffset
        jglo = moffset + this%dis%nodes + this%ioffset + n
        this%idxsymdglo(ipos) = matrix_sln%get_position_diag(iglo)
        this%idxsymoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
      end do
      !
      ! -- apt-apt contributions to gwf portion of global matrix
      if (this%idxbudfjf /= 0) then
        do ipos = 1, this%flowbudptr%budterm(this%idxbudfjf)%nlist
          n = this%flowbudptr%budterm(this%idxbudfjf)%id1(ipos)
          j = this%flowbudptr%budterm(this%idxbudfjf)%id2(ipos)
          iglo = moffset + this%dis%nodes + this%ioffset + n
          jglo = moffset + this%dis%nodes + this%ioffset + j
          this%idxfjfdglo(ipos) = matrix_sln%get_position_diag(iglo)
          this%idxfjfoffdglo(ipos) = matrix_sln%get_position(iglo, jglo)
        end do
      end if
    end if
  end subroutine apt_mc

  !> @brief Advanced package transport allocate and read (ar) routine
  !<
  subroutine apt_ar(this)
    ! -- modules
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    integer(I4B) :: j
    logical :: found
    ! -- formats
    character(len=*), parameter :: fmtapt = &
      "(1x,/1x,'APT -- ADVANCED PACKAGE TRANSPORT, VERSION 1, 3/5/2020', &
      &' INPUT READ FROM UNIT ', i0, //)"
    !
    ! -- Get obs setup
    call this%obs%obs_ar()
    !
    ! --print a message identifying the apt package.
    write (this%iout, fmtapt) this%inunit
    !
    ! -- Allocate arrays
    call this%apt_allocate_arrays()
    !
    ! -- read optional initial package parameters
    call this%read_initial_attr()
    !
    ! -- Find the package index in the GWF model or GWF budget file
    !    for the corresponding apt flow package
    call this%fmi%get_package_index(this%flowpackagename, this%igwfaptpak)
    !
    ! -- Tell fmi that this package is being handled by APT, otherwise
    !    SSM would handle the flows into GWT from this pack.  Then point the
    !    fmi data for an advanced package to xnewpak and qmfrommvr
    this%fmi%iatp(this%igwfaptpak) = 1
    this%fmi%datp(this%igwfaptpak)%concpack => this%get_mvr_depvar()
    this%fmi%datp(this%igwfaptpak)%qmfrommvr => this%qmfrommvr
    !
    ! -- If there is an associated flow package and the user wishes to put
    !    simulated concentrations (or temperatures) into a aux variable
    !    column, then find the column number.
    if (associated(this%flowpackagebnd)) then
      if (this%cauxfpconc /= '') then
        found = .false.
        do j = 1, this%flowpackagebnd%naux
          if (this%flowpackagebnd%auxname(j) == this%cauxfpconc) then
            this%iauxfpconc = j
            found = .true.
            exit
          end if
        end do
        if (this%iauxfpconc == 0) then
          errmsg = 'Could not find auxiliary variable '// &
                   trim(adjustl(this%cauxfpconc))//' in flow package '// &
                   trim(adjustl(this%flowpackagename))
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        else
          ! -- tell package not to update this auxiliary variable
          this%flowpackagebnd%noupdateauxvar(this%iauxfpconc) = 1
          call this%apt_copy2flowp()
        end if
      end if
    end if
  end subroutine apt_ar

  !> @brief Advanced package transport read and prepare (rp) routine
  !!
  !! This subroutine calls the attached packages' read and prepare routines.
  !<
  subroutine apt_rp(this)
    use TdisModule, only: kper, nper
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    integer(I4B) :: ierr
    integer(I4B) :: n
    logical :: isfound, endOfBlock
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: line
    integer(I4B) :: itemno
    integer(I4B) :: igwfnode
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Error.  Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    !
    ! -- set nbound to maxbound
    this%nbound = this%maxbound
    !
    ! -- Set ionper to the stress period number for which a new block of data
    !    will be read.
    if (this%inunit == 0) return
    !
    ! -- get stress period data
    if (this%ionper < kper) then
      !
      ! -- get period block
      call this%parser%GetBlock('PERIOD', isfound, ierr, &
                                supportOpenClose=.true., &
                                blockRequired=.false.)
      if (isfound) then
        !
        ! -- read ionper and check for increasing period numbers
        call this%read_check_ionper()
      else
        !
        ! -- PERIOD block not found
        if (ierr < 0) then
          ! -- End of file found; data applies for remainder of simulation.
          this%ionper = nper + 1
        else
          ! -- Found invalid block
          call this%parser%GetCurrentLine(line)
          write (errmsg, fmtblkerr) adjustl(trim(line))
          call store_error(errmsg)
          call this%parser%StoreErrorUnit()
        end if
      end if
    end if
    !
    ! -- Read data if ionper == kper
    if (this%ionper == kper) then
      !
      ! -- setup table for period data
      if (this%iprpak /= 0) then
        !
        ! -- reset the input table object
        title = trim(adjustl(this%text))//' PACKAGE ('// &
                trim(adjustl(this%packName))//') DATA FOR PERIOD'
        write (title, '(a,1x,i6)') trim(adjustl(title)), kper
        call table_cr(this%inputtab, this%packName, title)
        call this%inputtab%table_df(1, 4, this%iout, finalize=.FALSE.)
        text = 'NUMBER'
        call this%inputtab%initialize_column(text, 10, alignment=TABCENTER)
        text = 'KEYWORD'
        call this%inputtab%initialize_column(text, 20, alignment=TABLEFT)
        do n = 1, 2
          write (text, '(a,1x,i6)') 'VALUE', n
          call this%inputtab%initialize_column(text, 15, alignment=TABCENTER)
        end do
      end if
      !
      ! -- read data
      stressperiod: do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        !
        ! -- get feature number
        itemno = this%parser%GetInteger()
        !
        ! -- read data from the rest of the line
        call this%apt_set_stressperiod(itemno)
        !
        ! -- write line to table
        if (this%iprpak /= 0) then
          call this%parser%GetCurrentLine(line)
          call this%inputtab%line_to_columns(line)
        end if
      end do stressperiod

      if (this%iprpak /= 0) then
        call this%inputtab%finalize_table()
      end if
      !
      ! -- using stress period data from the previous stress period
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
    !
    ! -- write summary of stress period error messages
    ierr = count_errors()
    if (ierr > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- fill arrays
    do n = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(n)
      this%nodelist(n) = igwfnode
    end do
  end subroutine apt_rp

  subroutine apt_ad_chk(this)
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! function available for override by packages
  end subroutine apt_ad_chk

  !> @brief Advanced package transport set stress period routine.
  !!
  !! Set a stress period attribute for an advanced transport package feature
  !! (itemno) using keywords.
  !<
  subroutine apt_set_stressperiod(this, itemno)
    ! -- module
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    ! -- local
    character(len=LINELENGTH) :: text
    character(len=LINELENGTH) :: caux
    character(len=LINELENGTH) :: keyword
    integer(I4B) :: ierr
    integer(I4B) :: ii
    integer(I4B) :: jj
    real(DP), pointer :: bndElem => null()
    logical :: found
    ! -- formats
    !
    ! -- Support these general options in LKT, SFT, MWT, UZT
    ! STATUS <status>
    ! CONCENTRATION <concentration> or TEMPERATURE <temperature>
    ! WITHDRAWAL <withdrawal>
    ! AUXILIARY <auxname> <auxval>
    !
    ! -- read line
    call this%parser%GetStringCaps(keyword)
    select case (keyword)
    case ('STATUS')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetStringCaps(text)
      this%status(itemno) = text(1:8)
      if (text == 'CONSTANT') then
        this%iboundpak(itemno) = -1
      else if (text == 'INACTIVE') then
        this%iboundpak(itemno) = 0
      else if (text == 'ACTIVE') then
        this%iboundpak(itemno) = 1
      else
        write (errmsg, '(a,a)') &
          'Unknown '//trim(this%text)//' status keyword: ', text//'.'
        call store_error(errmsg)
      end if
    case ('CONCENTRATION', 'TEMPERATURE')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetString(text)
      jj = 1 ! For feature concentration
      bndElem => this%concfeat(itemno)
      call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                         this%packName, 'BND', this%tsManager, &
                                         this%iprpak, this%depvartype)
    case ('AUXILIARY')
      ierr = this%apt_check_valid(itemno)
      if (ierr /= 0) then
        goto 999
      end if
      call this%parser%GetStringCaps(caux)
      do jj = 1, this%naux
        if (trim(adjustl(caux)) /= trim(adjustl(this%auxname(jj)))) cycle
        call this%parser%GetString(text)
        ii = itemno
        bndElem => this%lauxvar(jj, ii)
        call read_value_or_time_series_adv(text, itemno, jj, bndElem, &
                                           this%packName, 'AUX', &
                                           this%tsManager, this%iprpak, &
                                           this%auxname(jj))
        exit
      end do
    case default
      !
      ! -- call the specific package to look for stress period data
      call this%pak_set_stressperiod(itemno, keyword, found)
      !
      ! -- terminate with error if data not valid
      if (.not. found) then
        write (errmsg, '(2a)') &
          'Unknown '//trim(adjustl(this%text))//' data keyword: ', &
          trim(keyword)//'.'
        call store_error(errmsg)
      end if
    end select
    !
    ! -- terminate if any errors were detected
999 if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine apt_set_stressperiod

  !> @brief Advanced package transport set stress period routine.
  !!
  !! Set a stress period attribute for an individual package. This routine
  !! must be overridden.
  !<
  subroutine pak_set_stressperiod(this, itemno, keyword, found)
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    character(len=*), intent(in) :: keyword
    logical, intent(inout) :: found
    ! -- local

    ! -- formats
    !
    ! -- this routine should never be called
    found = .false.
    call store_error('Program error: pak_set_stressperiod not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_set_stressperiod

  !> @brief Advanced package transport routine
  !!
  !! Determine if a valid feature number has been specified.
  !<
  function apt_check_valid(this, itemno) result(ierr)
    ! -- return
    integer(I4B) :: ierr
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    integer(I4B), intent(in) :: itemno
    ! -- formats
    ierr = 0
    if (itemno < 1 .or. itemno > this%ncv) then
      write (errmsg, '(a,1x,i6,1x,a,1x,i6)') &
        'Featureno ', itemno, 'must be > 0 and <= ', this%ncv
      call store_error(errmsg)
      ierr = 1
    end if
  end function apt_check_valid

  !> @brief Advanced package transport utility function
  !!
  !! Set the concentration (or temperature) to be used by either MVT or MVE
  !<
  function get_mvr_depvar(this)
    ! -- dummy
    class(TspAptType) :: this
    ! -- return
    real(dp), dimension(:), contiguous, pointer :: get_mvr_depvar
    !
    get_mvr_depvar => this%xnewpak
  end function get_mvr_depvar

  !> @brief Advanced package transport routine
  !!
  !! Add package connections to matrix
  !<
  subroutine apt_ad(this)
    ! -- modules
    use SimVariablesModule, only: iFailedStepRetry
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: j, iaux
    !
    ! -- Advance the time series
    call this%TsManager%ad()
    !
    ! -- update auxiliary variables by copying from the derived-type time
    !    series variable into the bndpackage auxvar variable so that this
    !    information is properly written to the GWF budget file
    if (this%naux > 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
        do iaux = 1, this%naux
          this%auxvar(iaux, j) = this%lauxvar(iaux, n)
        end do
      end do
    end if
    !
    ! -- copy xnew into xold and set xnewpak to specified concentration (or
    !    temperature) for constant concentration/temperature features
    if (iFailedStepRetry == 0) then
      do n = 1, this%ncv
        this%xoldpak(n) = this%xnewpak(n)
        if (this%iboundpak(n) < 0) then
          this%xnewpak(n) = this%concfeat(n)
        end if
      end do
    else
      do n = 1, this%ncv
        this%xnewpak(n) = this%xoldpak(n)
        if (this%iboundpak(n) < 0) then
          this%xnewpak(n) = this%concfeat(n)
        end if
      end do
    end if
    !
    ! -- pakmvrobj ad
    !if (this%imover == 1) then
    !  call this%pakmvrobj%ad()
    !end if
    !
    ! -- For each observation, push simulated value and corresponding
    !    simulation time from "current" to "preceding" and reset
    !    "current" value.
    call this%obs%obs_ad()
    !
    ! -- run package-specific checks
    call this%apt_ad_chk()
  end subroutine apt_ad

  !> @brief Override bnd reset for custom mover logic
  subroutine apt_reset(this)
    class(TspAptType) :: this !< GwtAptType object
    ! local
    integer(I4B) :: i
    !
    do i = 1, size(this%qmfrommvr)
      this%qmfrommvr(i) = DZERO
    end do
  end subroutine apt_reset

  subroutine apt_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    !
    ! -- Call fc depending on whether or not a matrix is expanded or not
    if (this%imatrows == 0) then
      call this%apt_fc_nonexpanded(rhs, ia, idxglo, matrix_sln)
    else
      call this%apt_fc_expanded(rhs, ia, idxglo, matrix_sln)
    end if
  end subroutine apt_fc

  !> @brief Advanced package transport fill coefficient (fc) method
  !!
  !! Routine to formulate the nonexpanded matrix case in which feature
  !! concentrations (or temperatures) are solved explicitly
  !<
  subroutine apt_fc_nonexpanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, igwfnode, idiag
    !
    ! -- solve for concentration (or temperatures) in the features
    call this%apt_solve()
    !
    ! -- add hcof and rhs terms (from apt_solve) to the gwf matrix
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(j)
      if (this%ibound(igwfnode) < 1) cycle
      idiag = idxglo(ia(igwfnode))
      call matrix_sln%add_value_pos(idiag, this%hcof(j))
      rhs(igwfnode) = rhs(igwfnode) + this%rhs(j)
    end do
  end subroutine apt_fc_nonexpanded

  !> @brief Advanced package transport fill coefficient (fc) method
  !!
  !! Routine to formulate the expanded matrix case in which new rows are added
  !! to the system of equations for each advanced package transport feature
  !<
  subroutine apt_fc_expanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: j, n, n1, n2
    integer(I4B) :: iloc
    integer(I4B) :: iposd, iposoffd
    integer(I4B) :: ipossymd, ipossymoffd
    real(DP) :: cold
    real(DP) :: qbnd, qbnd_scaled
    real(DP) :: omega
    real(DP) :: rrate
    real(DP) :: rhsval
    real(DP) :: hcofval
    !
    ! -- call the specific method for the advanced transport package, such as
    !    what would be overridden by
    !      GwtLktType, GwtSftType, GwtMwtType, GwtUztType
    !    This routine will add terms for rainfall, runoff, or other terms
    !    specific to the package
    call this%pak_fc_expanded(rhs, ia, idxglo, matrix_sln)
    !
    ! -- mass (or energy) storage in features
    do n = 1, this%ncv
      cold = this%xoldpak(n)
      iloc = this%idxlocnode(n)
      iposd = this%idxpakdiag(n)
      call this%apt_stor_term(n, n1, n2, rrate, rhsval, hcofval)
      call matrix_sln%add_value_pos(iposd, hcofval)
      rhs(iloc) = rhs(iloc) + rhsval
    end do
    !
    ! -- add to mover contribution
    if (this%idxbudtmvr /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudtmvr)%nlist
        call this%apt_tmvr_term(j, n1, n2, rrate, rhsval, hcofval)
        iloc = this%idxlocnode(n1)
        iposd = this%idxpakdiag(n1)
        call matrix_sln%add_value_pos(iposd, hcofval)
        rhs(iloc) = rhs(iloc) + rhsval
      end do
    end if
    !
    ! -- add from mover contribution
    if (this%idxbudfmvr /= 0) then
      do n = 1, this%ncv
        rhsval = this%qmfrommvr(n) ! this will already be in terms of energy for heat transport
        iloc = this%idxlocnode(n)
        rhs(iloc) = rhs(iloc) - rhsval
      end do
    end if
    !
    ! -- go through each apt-gwf connection
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      !
      ! -- set n to feature number and process if active feature
      n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      if (this%iboundpak(n) /= 0) then
        !
        ! -- set acoef and rhs to negative so they are relative to apt and not gwt
        qbnd = this%flowbudptr%budterm(this%idxbudgwf)%flow(j)
        omega = DZERO
        if (qbnd < DZERO) omega = DONE
        qbnd_scaled = qbnd * this%eqnsclfac
        !
        ! -- add to apt row
        iposd = this%idxdglo(j)
        iposoffd = this%idxoffdglo(j)
        call matrix_sln%add_value_pos(iposd, omega * qbnd_scaled)
        call matrix_sln%add_value_pos(iposoffd, (DONE - omega) * qbnd_scaled)
        !
        ! -- add to gwf row for apt connection
        ipossymd = this%idxsymdglo(j)
        ipossymoffd = this%idxsymoffdglo(j)
        call matrix_sln%add_value_pos(ipossymd, -(DONE - omega) * qbnd_scaled)
        call matrix_sln%add_value_pos(ipossymoffd, -omega * qbnd_scaled)
      end if
    end do
    !
    ! -- go through each apt-apt connection
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
        qbnd_scaled = qbnd * this%eqnsclfac
        iposd = this%idxfjfdglo(j)
        iposoffd = this%idxfjfoffdglo(j)
        call matrix_sln%add_value_pos(iposd, omega * qbnd_scaled)
        call matrix_sln%add_value_pos(iposoffd, (DONE - omega) * qbnd_scaled)
      end do
    end if
  end subroutine apt_fc_expanded

  !> @brief Advanced package transport fill coefficient (fc) method
  !!
  !! Routine to allow a subclass advanced transport package to inject
  !! terms into the matrix assembly.  This method must be overridden.
  !<
  subroutine pak_fc_expanded(this, rhs, ia, idxglo, matrix_sln)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_fc_expanded not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_fc_expanded

  !> @brief Advanced package transport routine
  !!
  !! Calculate advanced package transport hcof and rhs so transport budget is
  !! calculated.
  !<
  subroutine apt_cfupdate(this)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: j, n
    real(DP) :: qbnd
    real(DP) :: omega
    !
    ! -- Calculate hcof and rhs terms so GWF exchanges are calculated correctly
    ! -- go through each apt-gwf connection and calculate
    !    rhs and hcof terms for gwt/gwe matrix rows
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      this%hcof(j) = DZERO
      this%rhs(j) = DZERO
      if (this%iboundpak(n) /= 0) then
        qbnd = this%flowbudptr%budterm(this%idxbudgwf)%flow(j)
        omega = DZERO
        if (qbnd < DZERO) omega = DONE
        this%hcof(j) = -(DONE - omega) * qbnd * this%eqnsclfac
        this%rhs(j) = omega * qbnd * this%xnewpak(n) * this%eqnsclfac
      end if
    end do
  end subroutine apt_cfupdate

  !> @brief Advanced package transport calculate flows (cq) routine
  !!
  !! Calculate flows for the advanced package transport feature
  !<
  subroutine apt_cq(this, x, flowja, iadv)
    ! -- modules
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B), optional, intent(in) :: iadv
    ! -- local
    integer(I4B) :: n, n1, n2
    real(DP) :: rrate
    !
    ! -- Solve the feature concentrations (or temperatures) again or update
    !    the feature hcof and rhs terms
    if (this%imatrows == 0) then
      call this%apt_solve()
    else
      call this%apt_cfupdate()
    end if
    !
    ! -- call base functionality in bnd_cq
    call this%BndType%bnd_cq(x, flowja)
    !
    ! -- calculate storage term
    do n = 1, this%ncv
      rrate = DZERO
      if (this%iboundpak(n) > 0) then
        call this%apt_stor_term(n, n1, n2, rrate)
      end if
      this%qsto(n) = rrate
    end do
    !
    ! -- Copy concentrations (or temperatures) into the flow package auxiliary variable
    call this%apt_copy2flowp()
    !
    ! -- fill the budget object
    call this%apt_fill_budobj(x, flowja)
  end subroutine apt_cq

  !> @brief Save advanced package flows routine
  !<
  subroutine apt_ot_package_flows(this, icbcfl, ibudfl)
    use TdisModule, only: kstp, kper, delt, pertim, totim
    class(TspAptType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B) :: ibinun
    !
    ! -- write the flows from the budobj
    ibinun = 0
    if (this%ibudgetout /= 0) then
      ibinun = this%ibudgetout
    end if
    if (icbcfl == 0) ibinun = 0
    if (ibinun > 0) then
      call this%budobj%save_flows(this%dis, ibinun, kstp, kper, delt, &
                                  pertim, totim, this%iout)
    end if
    !
    ! -- Print lake flows table
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%budobj%write_flowtable(this%dis, kstp, kper)
    end if
  end subroutine apt_ot_package_flows

  subroutine apt_ot_dv(this, idvsave, idvprint)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    use TdisModule, only: kstp, kper, pertim, totim
    use ConstantsModule, only: DHNOFLO, DHDRY, LENBUDTXT
    use InputOutputModule, only: ulasav
    ! -- dummy
    class(TspAptType) :: this
    integer(I4B), intent(in) :: idvsave
    integer(I4B), intent(in) :: idvprint
    ! -- local
    integer(I4B) :: ibinun
    integer(I4B) :: n
    real(DP) :: c
    character(len=LENBUDTXT) :: text
    !
    ! -- set unit number for binary dependent variable output
    ibinun = 0
    if (this%iconcout /= 0) then
      ibinun = this%iconcout
    end if
    if (idvsave == 0) ibinun = 0
    !
    ! -- write binary output
    if (ibinun > 0) then
      do n = 1, this%ncv
        c = this%xnewpak(n)
        if (this%iboundpak(n) == 0) then
          c = DHNOFLO
        end if
        this%dbuff(n) = c
      end do
      write (text, '(a)') str_pad_left(this%depvartype, LENVARNAME)
      call ulasav(this%dbuff, text, kstp, kper, pertim, totim, &
                  this%ncv, 1, 1, ibinun)
    end if
    !
    ! -- write apt conc table
    if (idvprint /= 0 .and. this%iprconc /= 0) then
      !
      ! -- set table kstp and kper
      call this%dvtab%set_kstpkper(kstp, kper)
      !
      ! -- fill concentration data
      do n = 1, this%ncv
        if (this%inamedbound == 1) then
          call this%dvtab%add_term(this%featname(n))
        end if
        call this%dvtab%add_term(n)
        call this%dvtab%add_term(this%xnewpak(n))
      end do
    end if
  end subroutine apt_ot_dv

  !> @brief Print advanced package transport dependent variables
  !<
  subroutine apt_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- module
    use TdisModule, only: totim, delt
    ! -- dummy
    class(TspAptType) :: this !< TspAptType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    !
    call this%budobj%write_budtable(kstp, kper, iout, ibudfl, totim, delt)
  end subroutine apt_ot_bdsummary

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables for an advanced package
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    !
    ! -- allocate scalars in NumericalPackageType
    call this%BndType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%iauxfpconc, 'IAUXFPCONC', this%memoryPath)
    call mem_allocate(this%imatrows, 'IMATROWS', this%memoryPath)
    call mem_allocate(this%iprconc, 'IPRCONC', this%memoryPath)
    call mem_allocate(this%iconcout, 'ICONCOUT', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%igwfaptpak, 'IGWFAPTPAK', this%memoryPath)
    call mem_allocate(this%ncv, 'NCV', this%memoryPath)
    call mem_allocate(this%idxbudfjf, 'IDXBUDFJF', this%memoryPath)
    call mem_allocate(this%idxbudgwf, 'IDXBUDGWF', this%memoryPath)
    call mem_allocate(this%idxbudsto, 'IDXBUDSTO', this%memoryPath)
    call mem_allocate(this%idxbudtmvr, 'IDXBUDTMVR', this%memoryPath)
    call mem_allocate(this%idxbudfmvr, 'IDXBUDFMVR', this%memoryPath)
    call mem_allocate(this%idxbudaux, 'IDXBUDAUX', this%memoryPath)
    call mem_allocate(this%nconcbudssm, 'NCONCBUDSSM', this%memoryPath)
    call mem_allocate(this%idxprepak, 'IDXPREPAK', this%memoryPath)
    call mem_allocate(this%idxlastpak, 'IDXLASTPAK', this%memoryPath)
    !
    ! -- Initialize
    this%iauxfpconc = 0
    this%imatrows = 1
    this%iprconc = 0
    this%iconcout = 0
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%igwfaptpak = 0
    this%ncv = 0
    this%idxbudfjf = 0
    this%idxbudgwf = 0
    this%idxbudsto = 0
    this%idxbudtmvr = 0
    this%idxbudfmvr = 0
    this%idxbudaux = 0
    this%nconcbudssm = 0
    this%idxprepak = 0
    this%idxlastpak = 0
    !
    ! -- set this package as causing asymmetric matrix terms
    this%iasym = 1
  end subroutine allocate_scalars

  !> @ brief Allocate index arrays
  !!
  !! Allocate arrays that map to locations in the numerical solution
  !<
  subroutine apt_allocate_index_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    if (this%imatrows /= 0) then
      !
      ! -- count number of flow-ja-face connections
      n = 0
      if (this%idxbudfjf /= 0) then
        n = this%flowbudptr%budterm(this%idxbudfjf)%maxlist
      end if
      !
      ! -- allocate pointers to global matrix
      call mem_allocate(this%idxlocnode, this%ncv, 'IDXLOCNODE', &
                        this%memoryPath)
      call mem_allocate(this%idxpakdiag, this%ncv, 'IDXPAKDIAG', &
                        this%memoryPath)
      call mem_allocate(this%idxdglo, this%maxbound, 'IDXGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxoffdglo, this%maxbound, 'IDXOFFDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxsymdglo, this%maxbound, 'IDXSYMDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxsymoffdglo, this%maxbound, 'IDXSYMOFFDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxfjfdglo, n, 'IDXFJFDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxfjfoffdglo, n, 'IDXFJFOFFDGLO', &
                        this%memoryPath)
    else
      call mem_allocate(this%idxlocnode, 0, 'IDXLOCNODE', &
                        this%memoryPath)
      call mem_allocate(this%idxpakdiag, 0, 'IDXPAKDIAG', &
                        this%memoryPath)
      call mem_allocate(this%idxdglo, 0, 'IDXGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxoffdglo, 0, 'IDXOFFDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxsymdglo, 0, 'IDXSYMDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxsymoffdglo, 0, 'IDXSYMOFFDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxfjfdglo, 0, 'IDXFJFDGLO', &
                        this%memoryPath)
      call mem_allocate(this%idxfjfoffdglo, 0, 'IDXFJFOFFDGLO', &
                        this%memoryPath)
    end if
  end subroutine apt_allocate_index_arrays

  !> @ brief Allocate arrays
  !!
  !! Allocate advanced package transport arrays
  !<
  subroutine apt_allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- call standard BndType allocate scalars
    call this%BndType%allocate_arrays()
    !
    ! -- Allocate
    !
    ! -- allocate and initialize dbuff
    if (this%iconcout > 0) then
      call mem_allocate(this%dbuff, this%ncv, 'DBUFF', this%memoryPath)
      do n = 1, this%ncv
        this%dbuff(n) = DZERO
      end do
    else
      call mem_allocate(this%dbuff, 0, 'DBUFF', this%memoryPath)
    end if
    !
    ! -- allocate character array for status
    allocate (this%status(this%ncv))
    !
    ! -- time series
    call mem_allocate(this%concfeat, this%ncv, 'CONCFEAT', this%memoryPath)
    !
    ! -- budget terms
    call mem_allocate(this%qsto, this%ncv, 'QSTO', this%memoryPath)
    call mem_allocate(this%ccterm, this%ncv, 'CCTERM', this%memoryPath)
    !
    ! -- concentration for budget terms
    call mem_allocate(this%concbudssm, this%nconcbudssm, this%ncv, &
                      'CONCBUDSSM', this%memoryPath)
    !
    ! -- mass (or energy) added from the mover transport package
    call mem_allocate(this%qmfrommvr, this%ncv, 'QMFROMMVR', this%memoryPath)
    !
    ! -- initialize arrays
    do n = 1, this%ncv
      this%status(n) = 'ACTIVE'
      this%qsto(n) = DZERO
      this%ccterm(n) = DZERO
      this%qmfrommvr(n) = DZERO
      this%concbudssm(:, n) = DZERO
      this%concfeat(n) = DZERO
    end do
  end subroutine apt_allocate_arrays

  !> @ brief Deallocate memory
  !!
  !! Deallocate memory associated with this package
  !<
  subroutine apt_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    !
    ! -- deallocate arrays
    call mem_deallocate(this%dbuff)
    call mem_deallocate(this%qsto)
    call mem_deallocate(this%ccterm)
    call mem_deallocate(this%strt)
    call mem_deallocate(this%lauxvar)
    call mem_deallocate(this%xoldpak)
    if (this%imatrows == 0) then
      call mem_deallocate(this%iboundpak)
      call mem_deallocate(this%xnewpak)
    end if
    call mem_deallocate(this%concbudssm)
    call mem_deallocate(this%concfeat)
    call mem_deallocate(this%qmfrommvr)
    deallocate (this%status)
    deallocate (this%featname)
    !
    ! -- budobj
    call this%budobj%budgetobject_da()
    deallocate (this%budobj)
    nullify (this%budobj)
    !
    ! -- conc table
    if (this%iprconc > 0) then
      call this%dvtab%table_da()
      deallocate (this%dvtab)
      nullify (this%dvtab)
    end if
    !
    ! -- index pointers
    call mem_deallocate(this%idxlocnode)
    call mem_deallocate(this%idxpakdiag)
    call mem_deallocate(this%idxdglo)
    call mem_deallocate(this%idxoffdglo)
    call mem_deallocate(this%idxsymdglo)
    call mem_deallocate(this%idxsymoffdglo)
    call mem_deallocate(this%idxfjfdglo)
    call mem_deallocate(this%idxfjfoffdglo)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%iauxfpconc)
    call mem_deallocate(this%imatrows)
    call mem_deallocate(this%iprconc)
    call mem_deallocate(this%iconcout)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%igwfaptpak)
    call mem_deallocate(this%ncv)
    call mem_deallocate(this%idxbudfjf)
    call mem_deallocate(this%idxbudgwf)
    call mem_deallocate(this%idxbudsto)
    call mem_deallocate(this%idxbudtmvr)
    call mem_deallocate(this%idxbudfmvr)
    call mem_deallocate(this%idxbudaux)
    call mem_deallocate(this%idxbudssm)
    call mem_deallocate(this%nconcbudssm)
    call mem_deallocate(this%idxprepak)
    call mem_deallocate(this%idxlastpak)
    !
    ! -- deallocate scalars in NumericalPackageType
    call this%BndType%bnd_da()
  end subroutine apt_da

  !> @brief Find corresponding advanced package transport package
  !<
  subroutine find_apt_package(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_solve not implemented.', &
                     terminate=.TRUE.)
  end subroutine find_apt_package

  !> @brief Set options specific to the TspAptType
  !!
  !! This routine overrides BndType%bnd_options
  !<
  subroutine apt_options(this, option, found)
    use ConstantsModule, only: MAXCHARLEN, DZERO
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: urword, getunit, assign_iounit, openfile
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    character(len=*), intent(inout) :: option
    logical, intent(inout) :: found
    ! -- local
    character(len=MAXCHARLEN) :: fname, keyword
    ! -- formats
    character(len=*), parameter :: fmtaptbin = &
      "(4x, a, 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, &
      &/4x, 'OPENED ON UNIT: ', I0)"
    !
    found = .true.
    select case (option)
    case ('FLOW_PACKAGE_NAME')
      call this%parser%GetStringCaps(this%flowpackagename)
      write (this%iout, '(4x,a)') &
        'THIS '//trim(adjustl(this%text))//' PACKAGE CORRESPONDS TO A GWF &
        &PACKAGE WITH THE NAME '//trim(adjustl(this%flowpackagename))
    case ('FLOW_PACKAGE_AUXILIARY_NAME')
      call this%parser%GetStringCaps(this%cauxfpconc)
      write (this%iout, '(4x,a)') &
        'SIMULATED CONCENTRATIONS WILL BE COPIED INTO THE FLOW PACKAGE &
        &AUXILIARY VARIABLE WITH THE NAME '//trim(adjustl(this%cauxfpconc))
    case ('DEV_NONEXPANDING_MATRIX')
      ! -- use an iterative solution where concentration is not solved
      !    as part of the matrix.  It is instead solved separately with a
      !    general mixing equation and then added to the RHS of the GWT
      !    equations
      call this%parser%DevOpt()
      this%imatrows = 0
      write (this%iout, '(4x,a)') &
        trim(adjustl(this%text))// &
        ' WILL NOT ADD ADDITIONAL ROWS TO THE A MATRIX.'
    case ('PRINT_CONCENTRATION', 'PRINT_TEMPERATURE')
      this%iprconc = 1
      write (this%iout, '(4x,a,1x,a,1x,a)') trim(adjustl(this%text))// &
        trim(adjustl(this%depvartype))//'S WILL BE PRINTED TO LISTING &
          &FILE.'
    case ('CONCENTRATION', 'TEMPERATURE')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        this%iconcout = getunit()
        call openfile(this%iconcout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE')
        write (this%iout, fmtaptbin) &
          trim(adjustl(this%text)), trim(adjustl(this%depvartype)), &
          trim(fname), this%iconcout
      else
        write (errmsg, "('Optional', 1x, a, 1X, 'keyword must &
                         &be followed by FILEOUT')") this%depvartype
        call store_error(errmsg)
      end if
    case ('BUDGET')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
        call openfile(this%ibudgetout, this%iout, fname, 'DATA(BINARY)', &
                      form, access, 'REPLACE')
        write (this%iout, fmtaptbin) trim(adjustl(this%text)), 'BUDGET', &
          trim(fname), this%ibudgetout
      else
        call store_error('Optional BUDGET keyword must be followed by FILEOUT')
      end if
    case ('BUDGETCSV')
      call this%parser%GetStringCaps(keyword)
      if (keyword == 'FILEOUT') then
        call this%parser%GetString(fname)
        call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
        call openfile(this%ibudcsv, this%iout, fname, 'CSV', &
                      filstat_opt='REPLACE')
        write (this%iout, fmtaptbin) trim(adjustl(this%text)), 'BUDGET CSV', &
          trim(fname), this%ibudcsv
      else
        call store_error('Optional BUDGETCSV keyword must be followed by &
          &FILEOUT')
      end if
    case default
      !
      ! -- No options found
      found = .false.
    end select
  end subroutine apt_options

  !> @brief Determine dimensions for this advanced package
  !<
  subroutine apt_read_dimensions(this)
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    integer(I4B) :: ierr
    ! -- format
    !
    ! -- Set a pointer to the GWF LAK Package budobj
    if (this%flowpackagename == '') then
      this%flowpackagename = this%packName
      write (this%iout, '(4x,a)') &
        'THE FLOW PACKAGE NAME FOR '//trim(adjustl(this%text))//' WAS NOT &
        &SPECIFIED.  SETTING FLOW PACKAGE NAME TO '// &
        &trim(adjustl(this%flowpackagename))

    end if
    call this%find_apt_package()
    !
    ! -- Set dimensions from the GWF advanced package
    this%ncv = this%flowbudptr%ncv
    this%maxbound = this%flowbudptr%budterm(this%idxbudgwf)%maxlist
    this%nbound = this%maxbound
    write (this%iout, '(a, a)') 'SETTING DIMENSIONS FOR PACKAGE ', this%packName
    write (this%iout, '(2x,a,i0)') 'NUMBER OF CONTROL VOLUMES = ', this%ncv
    write (this%iout, '(2x,a,i0)') 'MAXBOUND = ', this%maxbound
    write (this%iout, '(2x,a,i0)') 'NBOUND = ', this%nbound
    if (this%imatrows /= 0) then
      this%npakeq = this%ncv
      write (this%iout, '(2x,a)') trim(adjustl(this%text))// &
        ' SOLVED AS PART OF GWT MATRIX EQUATIONS'
    else
      write (this%iout, '(2x,a)') trim(adjustl(this%text))// &
        ' SOLVED SEPARATELY FROM GWT MATRIX EQUATIONS '
    end if
    write (this%iout, '(a, //)') 'DONE SETTING DIMENSIONS FOR '// &
      trim(adjustl(this%text))
    !
    ! -- Check for errors
    if (this%ncv < 0) then
      write (errmsg, '(a)') &
        'Number of control volumes could not be determined correctly.'
      call store_error(errmsg)
    end if
    !
    ! -- stop if errors were encountered in the DIMENSIONS block
    ierr = count_errors()
    if (ierr > 0) then
      call store_error_unit(this%inunit)
    end if
    !
    ! -- read packagedata block
    call this%apt_read_cvs()
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
    !
    ! -- setup the budget object
    call this%apt_setup_budobj()
    !
    ! -- setup the conc table object
    call this%apt_setup_tableobj()
  end subroutine apt_read_dimensions

  !> @brief Read feature information for this advanced package
  !<
  subroutine apt_read_cvs(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use TimeSeriesManagerModule, only: read_value_or_time_series_adv
    ! -- dummy
    class(TspAptType), intent(inout) :: this
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
  end subroutine apt_read_cvs

  !> @brief Read the initial parameters for an advanced package
  !<
  subroutine apt_read_initial_attr(this)
    use ConstantsModule, only: LINELENGTH
    use BudgetModule, only: budget_cr
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    !character(len=LINELENGTH) :: text
    integer(I4B) :: j, n

    !
    ! -- initialize xnewpak and set feature concentration (or temperature)
    ! -- todo: this should be a time series?
    do n = 1, this%ncv
      this%xnewpak(n) = this%strt(n)
      !
      ! -- todo: read aux
      !
      ! -- todo: read boundname
    end do
    !
    ! -- initialize status (iboundpak) of lakes to active
    do n = 1, this%ncv
      if (this%status(n) == 'CONSTANT') then
        this%iboundpak(n) = -1
      else if (this%status(n) == 'INACTIVE') then
        this%iboundpak(n) = 0
      else if (this%status(n) == 'ACTIVE ') then
        this%iboundpak(n) = 1
      end if
    end do
    !
    ! -- set boundname for each connection
    if (this%inamedbound /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
        this%boundname(j) = this%featname(n)
      end do
    end if
    !
    ! -- copy boundname into boundname_cst
    call this%copy_boundname()
  end subroutine apt_read_initial_attr

  !> @brief Add terms specific to advanced package transport to the explicit
  !! solve
  !!
  !! Explicit solve for concentration (or temperature) in advaced package
  !! features, which is an alternative to the iterative implicit solve.
  !<
  subroutine apt_solve(this)
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: n, j, igwfnode
    integer(I4B) :: n1, n2
    real(DP) :: rrate
    real(DP) :: ctmp
    real(DP) :: c1, qbnd
    real(DP) :: hcofval, rhsval
    !
    ! -- initialize dbuff
    do n = 1, this%ncv
      this%dbuff(n) = DZERO
    end do
    !
    ! -- call the individual package routines to add terms specific to the
    !    advanced transport package
    call this%pak_solve()
    !
    ! -- add to mover contribution
    if (this%idxbudtmvr /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudtmvr)%nlist
        call this%apt_tmvr_term(j, n1, n2, rrate)
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- add from mover contribution
    if (this%idxbudfmvr /= 0) then
      do n1 = 1, size(this%qmfrommvr)
        rrate = this%qmfrommvr(n1) ! Will be in terms of energy for heat transport
        this%dbuff(n1) = this%dbuff(n1) + rrate
      end do
    end if
    !
    ! -- go through each gwf connection and accumulate
    !    total mass (or energy) in dbuff mass
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      this%hcof(j) = DZERO
      this%rhs(j) = DZERO
      igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(j)
      qbnd = this%flowbudptr%budterm(this%idxbudgwf)%flow(j)
      if (qbnd <= DZERO) then
        ctmp = this%xnewpak(n)
        this%rhs(j) = qbnd * ctmp * this%eqnsclfac
      else
        ctmp = this%xnew(igwfnode)
        this%hcof(j) = -qbnd * this%eqnsclfac
      end if
      c1 = qbnd * ctmp * this%eqnsclfac
      this%dbuff(n) = this%dbuff(n) + c1
    end do
    !
    ! -- go through each "within apt-apt" connection (e.g., lak-lak) and
    !    accumulate total mass (or energy) in dbuff mass
    if (this%idxbudfjf /= 0) then
      do j = 1, this%flowbudptr%budterm(this%idxbudfjf)%nlist
        call this%apt_fjf_term(j, n1, n2, rrate)
        c1 = rrate
        this%dbuff(n1) = this%dbuff(n1) + c1
      end do
    end if
    !
    ! -- calculate the feature concentration/temperature
    do n = 1, this%ncv
      call this%apt_stor_term(n, n1, n2, rrate, rhsval, hcofval)
      !
      ! -- at this point, dbuff has q * c for all sources, so now
      !    add Vold / dt * Cold
      this%dbuff(n) = this%dbuff(n) - rhsval
      !
      ! -- Now to calculate c, need to divide dbuff by hcofval
      c1 = -this%dbuff(n) / hcofval
      if (this%iboundpak(n) > 0) then
        this%xnewpak(n) = c1
      end if
    end do
  end subroutine apt_solve

  !> @brief Add terms specific to advanced package transport features to the
  !! explicit solve routine
  !!
  !! This routine must be overridden by the specific apt package
  !<
  subroutine pak_solve(this)
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_solve not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_solve

  !> @brief Accumulate constant concentration (or temperature) terms for budget
  !<
  subroutine apt_accumulate_ccterm(this, ilak, rrate, ccratin, ccratout)
    ! -- dummy
    class(TspAptType) :: this
    integer(I4B), intent(in) :: ilak
    real(DP), intent(in) :: rrate
    real(DP), intent(inout) :: ccratin
    real(DP), intent(inout) :: ccratout
    ! -- locals
    real(DP) :: q
    ! format
    ! code
    !
    if (this%iboundpak(ilak) < 0) then
      q = -rrate
      this%ccterm(ilak) = this%ccterm(ilak) + q
      !
      ! -- See if flow is into lake or out of lake.
      if (q < DZERO) then
        !
        ! -- Flow is out of lake subtract rate from ratout.
        ccratout = ccratout - q
      else
        !
        ! -- Flow is into lake; add rate to ratin.
        ccratin = ccratin + q
      end if
    end if
  end subroutine apt_accumulate_ccterm

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used.
  !<
  subroutine define_listlabel(this)
    class(TspAptType), intent(inout) :: this
    !
    ! -- create the header list label
    this%listlabel = trim(this%filtyp)//' NO.'
    if (this%dis%ndim == 3) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'ROW'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'COL'
    elseif (this%dis%ndim == 2) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'CELL2D'
    else
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'NODE'
    end if
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'STRESS RATE'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine define_listlabel

  !> @brief Set pointers to model arrays and variables so that a package has
  !! access to these items.
  !<
  subroutine apt_set_pointers(this, neq, ibound, xnew, xold, flowja)
    class(TspAptType) :: this
    integer(I4B), pointer :: neq
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    real(DP), dimension(:), pointer, contiguous :: xnew
    real(DP), dimension(:), pointer, contiguous :: xold
    real(DP), dimension(:), pointer, contiguous :: flowja
    ! -- local
    integer(I4B) :: istart, iend
    !
    ! -- call base BndType set_pointers
    call this%BndType%set_pointers(neq, ibound, xnew, xold, flowja)
    !
    ! -- Set the pointers
    !
    ! -- set package pointers
    if (this%imatrows /= 0) then
      istart = this%dis%nodes + this%ioffset + 1
      iend = istart + this%ncv - 1
      this%iboundpak => this%ibound(istart:iend)
      this%xnewpak => this%xnew(istart:iend)
    end if
  end subroutine apt_set_pointers

  !> @brief Return the feature new volume and old volume
  !<
  subroutine apt_get_volumes(this, icv, vnew, vold, delt)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    integer(I4B), intent(in) :: icv
    real(DP), intent(inout) :: vnew, vold
    real(DP), intent(in) :: delt
    ! -- local
    real(DP) :: qss
    !
    ! -- get volumes
    vold = DZERO
    vnew = vold
    if (this%idxbudsto /= 0) then
      qss = this%flowbudptr%budterm(this%idxbudsto)%flow(icv)
      vnew = this%flowbudptr%budterm(this%idxbudsto)%auxvar(1, icv)
      vold = vnew + qss * delt
    end if
  end subroutine apt_get_volumes

  !> @brief Function to return the number of budget terms just for this package
  !!
  !! This function must be overridden.
  !<
  function pak_get_nbudterms(this) result(nbudterms)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    ! -- return
    integer(I4B) :: nbudterms
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_get_nbudterms not implemented.', &
                     terminate=.TRUE.)
    nbudterms = 0
  end function pak_get_nbudterms

  !> @brief Set up the budget object that stores advanced package flow terms
  !<
  subroutine apt_setup_budobj(this)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: nlen
    integer(I4B) :: n, n1, n2
    integer(I4B) :: maxlist, naux
    integer(I4B) :: idx
    logical :: ordered_id1
    real(DP) :: q
    character(len=LENBUDTXT) :: bddim_opt
    character(len=LENBUDTXT) :: text, textt
    character(len=LENBUDTXT), dimension(1) :: auxtxt
    !
    ! -- initialize nbudterm
    nbudterm = 0
    !
    ! -- Determine if there are flow-ja-face terms
    nlen = 0
    if (this%idxbudfjf /= 0) then
      nlen = this%flowbudptr%budterm(this%idxbudfjf)%maxlist
    end if
    !
    ! -- Determine the number of budget terms associated with apt.
    !    These are fixed for the simulation and cannot change
    !
    ! -- add one if flow-ja-face present
    if (this%idxbudfjf /= 0) nbudterm = nbudterm + 1
    !
    ! -- All the APT packages have GWF, STORAGE, and CONSTANT
    nbudterm = nbudterm + 3
    !
    ! -- add terms for the specific package
    nbudterm = nbudterm + this%pak_get_nbudterms()
    !
    ! -- add for mover terms and auxiliary
    if (this%idxbudtmvr /= 0) nbudterm = nbudterm + 1
    if (this%idxbudfmvr /= 0) nbudterm = nbudterm + 1
    if (this%naux > 0) nbudterm = nbudterm + 1
    !
    ! -- set up budobj
    call budgetobject_cr(this%budobj, this%packName)
    !
    bddim_opt = this%depvarunitabbrev
    call this%budobj%budgetobject_df(this%ncv, nbudterm, 0, 0, &
                                     bddim_opt=bddim_opt, ibudcsv=this%ibudcsv)
    idx = 0
    !
    ! -- Go through and set up each budget term
    if (nlen > 0) then
      text = '    FLOW-JA-FACE'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudfjf)%maxlist
      naux = 0
      ordered_id1 = this%flowbudptr%budterm(this%idxbudfjf)%ordered_id1
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, ordered_id1=ordered_id1)
      !
      ! -- store outlet connectivity
      call this%budobj%budterm(idx)%reset(maxlist)
      q = DZERO
      do n = 1, maxlist
        n1 = this%flowbudptr%budterm(this%idxbudfjf)%id1(n)
        n2 = this%flowbudptr%budterm(this%idxbudfjf)%id2(n)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
      end do
    end if
    !
    ! --
    text = '             GWF'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudgwf)%maxlist
    naux = 0
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%name_model, &
                                             maxlist, .false., .true., &
                                             naux)
    call this%budobj%budterm(idx)%reset(maxlist)
    q = DZERO
    do n = 1, maxlist
      n1 = this%flowbudptr%budterm(this%idxbudgwf)%id1(n)
      n2 = this%flowbudptr%budterm(this%idxbudgwf)%id2(n)
      call this%budobj%budterm(idx)%update_term(n1, n2, q)
    end do
    !
    ! -- Reserve space for the package specific terms
    this%idxprepak = idx
    call this%pak_setup_budobj(idx)
    this%idxlastpak = idx
    !
    ! --
    text = '         STORAGE'
    idx = idx + 1
    maxlist = this%flowbudptr%budterm(this%idxbudsto)%maxlist
    naux = 1
    write (textt, '(a)') str_pad_left(this%depvarunit, 16)
    auxtxt(1) = textt ! '            MASS' or '          ENERGY'
    call this%budobj%budterm(idx)%initialize(text, &
                                             this%name_model, &
                                             this%packName, &
                                             this%name_model, &
                                             this%packName, &
                                             maxlist, .false., .false., &
                                             naux, auxtxt)
    if (this%idxbudtmvr /= 0) then
      !
      ! --
      text = '          TO-MVR'
      idx = idx + 1
      maxlist = this%flowbudptr%budterm(this%idxbudtmvr)%maxlist
      naux = 0
      ordered_id1 = this%flowbudptr%budterm(this%idxbudtmvr)%ordered_id1
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, ordered_id1=ordered_id1)
    end if
    if (this%idxbudfmvr /= 0) then
      !
      ! --
      text = '        FROM-MVR'
      idx = idx + 1
      maxlist = this%ncv
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
    text = '        CONSTANT'
    idx = idx + 1
    maxlist = this%ncv
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
    naux = this%naux
    if (naux > 0) then
      !
      ! --
      text = '       AUXILIARY'
      idx = idx + 1
      maxlist = this%ncv
      call this%budobj%budterm(idx)%initialize(text, &
                                               this%name_model, &
                                               this%packName, &
                                               this%name_model, &
                                               this%packName, &
                                               maxlist, .false., .false., &
                                               naux, this%auxname)
    end if
    !
    ! -- if flow for each control volume are written to the listing file
    if (this%iprflow /= 0) then
      call this%budobj%flowtable_df(this%iout)
    end if
  end subroutine apt_setup_budobj

  !> @brief Set up a budget object that stores an advanced package flows
  !!
  !! Individual packages set up their budget terms.  Must be overridden.
  !<
  subroutine pak_setup_budobj(this, idx)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    integer(I4B), intent(inout) :: idx
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_setup_budobj not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_setup_budobj

  !> @brief Copy flow terms into this%budobj
  !<
  subroutine apt_fill_budobj(this, x, flowja)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy
    class(TspAptType) :: this
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    ! -- local
    integer(I4B) :: naux
    real(DP), dimension(:), allocatable :: auxvartmp
    integer(I4B) :: i, j, n1, n2
    integer(I4B) :: idx
    integer(I4B) :: nlen
    integer(I4B) :: nlist
    integer(I4B) :: igwfnode
    real(DP) :: q
    real(DP) :: v0, v1
    real(DP) :: ccratin, ccratout
    ! -- formats
    !
    ! -- initialize counter
    idx = 0
    !
    ! -- initialize ccterm, which is used to sum up all mass (or energy) flows
    !    into a constant concentration (or temperature) cell
    ccratin = DZERO
    ccratout = DZERO
    do n1 = 1, this%ncv
      this%ccterm(n1) = DZERO
    end do
    !
    ! -- FLOW JA FACE
    nlen = 0
    if (this%idxbudfjf /= 0) then
      nlen = this%flowbudptr%budterm(this%idxbudfjf)%maxlist
    end if
    if (nlen > 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudfjf)%maxlist
      call this%budobj%budterm(idx)%reset(nlist)
      q = DZERO
      do j = 1, nlist
        call this%apt_fjf_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- GWF (LEAKAGE)
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%maxbound)
    do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
      q = DZERO
      n1 = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
      if (this%iboundpak(n1) /= 0) then
        igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(j)
        q = this%hcof(j) * x(igwfnode) - this%rhs(j)
        q = -q ! flip sign so relative to advanced package feature
      end if
      call this%budobj%budterm(idx)%update_term(n1, igwfnode, q)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    !
    ! -- skip individual package terms for now and process them last
    ! -- in case they depend on the other terms (as for uze)
    idx = this%idxlastpak
    !
    ! -- STORAGE
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%ncv)
    allocate (auxvartmp(1))
    do n1 = 1, this%ncv
      call this%apt_get_volumes(n1, v1, v0, delt)
      auxvartmp(1) = v1 * this%xnewpak(n1) ! Note: When GWE is added, check if this needs a factor of eqnsclfac
      q = this%qsto(n1)
      call this%budobj%budterm(idx)%update_term(n1, n1, q, auxvartmp)
      call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
    end do
    deallocate (auxvartmp)
    !
    ! -- TO MOVER
    if (this%idxbudtmvr /= 0) then
      idx = idx + 1
      nlist = this%flowbudptr%budterm(this%idxbudtmvr)%nlist
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%apt_tmvr_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n2, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- FROM MOVER
    if (this%idxbudfmvr /= 0) then
      idx = idx + 1
      nlist = this%ncv
      call this%budobj%budterm(idx)%reset(nlist)
      do j = 1, nlist
        call this%apt_fmvr_term(j, n1, n2, q)
        call this%budobj%budterm(idx)%update_term(n1, n1, q)
        call this%apt_accumulate_ccterm(n1, q, ccratin, ccratout)
      end do
    end if
    !
    ! -- CONSTANT FLOW
    idx = idx + 1
    call this%budobj%budterm(idx)%reset(this%ncv)
    do n1 = 1, this%ncv
      q = this%ccterm(n1)
      call this%budobj%budterm(idx)%update_term(n1, n1, q)
    end do
    !
    ! -- AUXILIARY VARIABLES
    naux = this%naux
    if (naux > 0) then
      idx = idx + 1
      allocate (auxvartmp(naux))
      call this%budobj%budterm(idx)%reset(this%ncv)
      do n1 = 1, this%ncv
        q = DZERO
        do i = 1, naux
          auxvartmp(i) = this%lauxvar(i, n1)
        end do
        call this%budobj%budterm(idx)%update_term(n1, n1, q, auxvartmp)
      end do
      deallocate (auxvartmp)
    end if
    !
    ! -- individual package terms processed last
    idx = this%idxprepak
    call this%pak_fill_budobj(idx, x, flowja, ccratin, ccratout)
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine apt_fill_budobj

  !> @brief Copy flow terms into this%budobj, must be overridden
  !<
  subroutine pak_fill_budobj(this, idx, x, flowja, ccratin, ccratout)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    integer(I4B), intent(inout) :: idx
    real(DP), dimension(:), intent(in) :: x
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    real(DP), intent(inout) :: ccratin
    real(DP), intent(inout) :: ccratout
    ! -- local
    ! -- formats
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_fill_budobj not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_fill_budobj

  !> @brief Account for mass or energy storage in advanced package features
  !<
  subroutine apt_stor_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    use TdisModule, only: delt
    class(TspAptType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    real(DP) :: v0, v1
    real(DP) :: c0, c1
    !
    n1 = ientry
    n2 = ientry
    call this%apt_get_volumes(n1, v1, v0, delt)
    c0 = this%xoldpak(n1)
    c1 = this%xnewpak(n1)
    if (present(rrate)) then
      rrate = (-c1 * v1 / delt + c0 * v0 / delt) * this%eqnsclfac
    end if
    !
    if (present(rhsval)) rhsval = -c0 * v0 * this%eqnsclfac / delt
    if (present(hcofval)) hcofval = -v1 * this%eqnsclfac / delt
  end subroutine apt_stor_term

  !> @brief Account for mass or energy transferred to the MVR package
  !<
  subroutine apt_tmvr_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
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
    ! -- Calculate MVR-related terms
    n1 = this%flowbudptr%budterm(this%idxbudtmvr)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudtmvr)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudtmvr)%flow(ientry)
    ctmp = this%xnewpak(n1)
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = DZERO
    if (present(hcofval)) hcofval = qbnd * this%eqnsclfac
  end subroutine apt_tmvr_term

  !> @brief Account for mass or energy transferred to this package from the
  !! MVR package
  !<
  subroutine apt_fmvr_term(this, ientry, n1, n2, rrate, &
                           rhsval, hcofval)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    integer(I4B), intent(in) :: ientry
    integer(I4B), intent(inout) :: n1
    integer(I4B), intent(inout) :: n2
    real(DP), intent(inout), optional :: rrate
    real(DP), intent(inout), optional :: rhsval
    real(DP), intent(inout), optional :: hcofval
    !
    ! -- Calculate MVR-related terms
    n1 = ientry
    n2 = n1
    if (present(rrate)) rrate = this%qmfrommvr(n1) ! NOTE: When bringing in GWE, ensure this is in terms of energy. Might need to apply eqnsclfac here.
    if (present(rhsval)) rhsval = this%qmfrommvr(n1)
    if (present(hcofval)) hcofval = DZERO
  end subroutine apt_fmvr_term

  !> @brief Go through each "within apt-apt" connection (e.g., lkt-lkt, or
  !! sft-sft) and accumulate total mass (or energy) in dbuff mass
  !<
  subroutine apt_fjf_term(this, ientry, n1, n2, rrate, &
                          rhsval, hcofval)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
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
    n1 = this%flowbudptr%budterm(this%idxbudfjf)%id1(ientry)
    n2 = this%flowbudptr%budterm(this%idxbudfjf)%id2(ientry)
    qbnd = this%flowbudptr%budterm(this%idxbudfjf)%flow(ientry)
    if (qbnd <= 0) then
      ctmp = this%xnewpak(n1)
    else
      ctmp = this%xnewpak(n2)
    end if
    if (present(rrate)) rrate = ctmp * qbnd * this%eqnsclfac
    if (present(rhsval)) rhsval = -rrate * this%eqnsclfac
    if (present(hcofval)) hcofval = DZERO
  end subroutine apt_fjf_term

  !> @brief Copy concentrations (or temperatures) into flow package aux
  !! variable
  !<
  subroutine apt_copy2flowp(this)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: n, j
    !
    ! -- copy
    if (this%iauxfpconc /= 0) then
      !
      ! -- go through each apt-gwf connection
      do j = 1, this%flowbudptr%budterm(this%idxbudgwf)%nlist
        !
        ! -- set n to feature number and process if active feature
        n = this%flowbudptr%budterm(this%idxbudgwf)%id1(j)
        this%flowpackagebnd%auxvar(this%iauxfpconc, j) = this%xnewpak(n)
      end do
    end if
  end subroutine apt_copy2flowp

  !> @brief Determine whether an obs type is supported
  !!
  !! This function:
  !!   - returns true if APT package supports named observation.
  !!   - overrides BndType%bnd_obs_supported()
  !<
  logical function apt_obs_supported(this)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    !
    ! -- Set to true
    apt_obs_supported = .true.
  end function apt_obs_supported

  !> @brief Define observation type
  !!
  !! This routine:
  !!   - stores observation types supported by APT package.
  !!   - overrides BndType%bnd_df_obs
  !<
  subroutine apt_df_obs(this)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    !
    ! -- call additional specific observations for lkt, sft, mwt, and uzt
    call this%pak_df_obs()
  end subroutine apt_df_obs

  !> @brief Define apt observation type
  !!
  !! This routine:
  !!   - stores observations supported by the APT package
  !!   - must be overridden by child class
  subroutine pak_df_obs(this)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_df_obs not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_df_obs

  !> @brief Process package specific obs
  !!
  !! Method to process specific observations for this package.
  !<
  subroutine pak_rp_obs(this, obsrv, found)
    ! -- dummy
    class(TspAptType), intent(inout) :: this !< package class
    type(ObserveType), intent(inout) :: obsrv !< observation object
    logical, intent(inout) :: found !< indicate whether observation was found
    ! -- local
    !
    ! -- this routine should never be called
    call store_error('Program error: pak_rp_obs not implemented.', &
                     terminate=.TRUE.)
  end subroutine pak_rp_obs

  !> @brief Prepare observation
  !!
  !! Find the indices for this observation assuming they are indexed by
  !! feature number
  !<
  subroutine rp_obs_byfeature(this, obsrv)
    class(TspAptType), intent(inout) :: this !< object
    type(ObserveType), intent(inout) :: obsrv !< observation
    integer(I4B) :: nn1
    integer(I4B) :: j
    logical :: jfound
    character(len=*), parameter :: fmterr = &
      "('Boundary ', a, ' for observation ', a, &
      &' is invalid in package ', a)"
    nn1 = obsrv%NodeNumber
    if (nn1 == NAMEDBOUNDFLAG) then
      jfound = .false.
      do j = 1, this%ncv
        if (this%featname(j) == obsrv%FeatureName) then
          jfound = .true.
          call obsrv%AddObsIndex(j)
        end if
      end do
      if (.not. jfound) then
        write (errmsg, fmterr) trim(obsrv%FeatureName), trim(obsrv%Name), &
          trim(this%packName)
        call store_error(errmsg)
      end if
    else
      !
      ! -- ensure nn1 is > 0 and < ncv
      if (nn1 < 0 .or. nn1 > this%ncv) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), ' was assigned ID = ', nn1, &
          '.  ID must be >= 1 and <= ', this%ncv, '.'
        call store_error(errmsg)
      end if
      call obsrv%AddObsIndex(nn1)
    end if
  end subroutine rp_obs_byfeature

  !> @brief Prepare observation
  !!
  !! Find the indices for this observation assuming they are first indexed
  !! by feature number and secondly by a connection number
  !<
  subroutine rp_obs_budterm(this, obsrv, budterm)
    class(TspAptType), intent(inout) :: this !< object
    type(ObserveType), intent(inout) :: obsrv !< observation
    type(BudgetTermType), intent(in) :: budterm !< budget term
    integer(I4B) :: nn1
    integer(I4B) :: iconn
    integer(I4B) :: icv
    integer(I4B) :: idx
    integer(I4B) :: j
    logical :: jfound
    character(len=*), parameter :: fmterr = &
      "('Boundary ', a, ' for observation ', a, &
      &' is invalid in package ', a)"
    nn1 = obsrv%NodeNumber
    if (nn1 == NAMEDBOUNDFLAG) then
      jfound = .false.
      do j = 1, budterm%nlist
        icv = budterm%id1(j)
        if (this%featname(icv) == obsrv%FeatureName) then
          jfound = .true.
          call obsrv%AddObsIndex(j)
        end if
      end do
      if (.not. jfound) then
        write (errmsg, fmterr) trim(obsrv%FeatureName), trim(obsrv%Name), &
          trim(this%packName)
        call store_error(errmsg)
      end if
    else
      !
      ! -- ensure nn1 is > 0 and < ncv
      if (nn1 < 0 .or. nn1 > this%ncv) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), ' was assigned ID = ', nn1, &
          '.  ID must be >= 1 and <= ', this%ncv, '.'
        call store_error(errmsg)
      end if
      iconn = obsrv%NodeNumber2
      do j = 1, budterm%nlist
        if (budterm%id1(j) == nn1) then
          ! -- Look for the first occurrence of nn1, then set indxbnds
          !    to the iconn record after that
          idx = j + iconn - 1
          call obsrv%AddObsIndex(idx)
          exit
        end if
      end do
      if (idx < 1 .or. idx > budterm%nlist) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), ' specifies iconn = ', iconn, &
          ',  but this is not a valid connection for ID ', nn1, '.'
        call store_error(errmsg)
      else if (budterm%id1(idx) /= nn1) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), ' specifies iconn = ', iconn, &
          ',  but this is not a valid connection for ID ', nn1, '.'
        call store_error(errmsg)
      end if
    end if
  end subroutine rp_obs_budterm

  !> @brief Prepare observation
  !!
  !! Find the indices for this observation assuming they are first indexed
  !! by a feature number and secondly by a second feature number
  !<
  subroutine rp_obs_flowjaface(this, obsrv, budterm)
    class(TspAptType), intent(inout) :: this !< object
    type(ObserveType), intent(inout) :: obsrv !< observation
    type(BudgetTermType), intent(in) :: budterm !< budget term
    integer(I4B) :: nn1
    integer(I4B) :: nn2
    integer(I4B) :: icv
    integer(I4B) :: j
    logical :: jfound
    character(len=*), parameter :: fmterr = &
      "('Boundary ', a, ' for observation ', a, &
      &' is invalid in package ', a)"
    nn1 = obsrv%NodeNumber
    if (nn1 == NAMEDBOUNDFLAG) then
      jfound = .false.
      do j = 1, budterm%nlist
        icv = budterm%id1(j)
        if (this%featname(icv) == obsrv%FeatureName) then
          jfound = .true.
          call obsrv%AddObsIndex(j)
        end if
      end do
      if (.not. jfound) then
        write (errmsg, fmterr) trim(obsrv%FeatureName), trim(obsrv%Name), &
          trim(this%packName)
        call store_error(errmsg)
      end if
    else
      !
      ! -- ensure nn1 is > 0 and < ncv
      if (nn1 < 0 .or. nn1 > this%ncv) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), ' was assigned ID = ', nn1, &
          '.  ID must be >= 1 and <= ', this%ncv, '.'
        call store_error(errmsg)
      end if
      nn2 = obsrv%NodeNumber2
      !
      ! -- ensure nn2 is > 0 and < ncv
      if (nn2 < 0 .or. nn2 > this%ncv) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), ' was assigned ID2 = ', nn2, &
          '.  ID must be >= 1 and <= ', this%ncv, '.'
        call store_error(errmsg)
      end if
      ! -- Look for nn1 and nn2 in id1 and id2
      jfound = .false.
      do j = 1, budterm%nlist
        if (budterm%id1(j) == nn1 .and. budterm%id2(j) == nn2) then
          call obsrv%AddObsIndex(j)
          jfound = .true.
        end if
      end do
      if (.not. jfound) then
        write (errmsg, '(7a, i0, a, i0, a)') &
          'Observation ', trim(obsrv%Name), ' of type ', &
          trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
          trim(this%packName), &
          ' specifies a connection between feature ', nn1, &
          ' feature ', nn2, ', but these features are not connected.'
        call store_error(errmsg)
      end if
    end if
  end subroutine rp_obs_flowjaface

  !> @brief Read and prepare apt-related observations
  !!
  !! Method to process specific observations for an apt package
  !<
  subroutine apt_rp_obs(this)
    ! -- modules
    use TdisModule, only: kper
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    logical :: found
    class(ObserveType), pointer :: obsrv => null()
    !
    if (kper == 1) then
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        select case (obsrv%ObsTypeId)
        case ('CONCENTRATION', 'TEMPERATURE')
          call this%rp_obs_byfeature(obsrv)
          !
          ! -- catch non-cumulative observation assigned to observation defined
          !    by a boundname that is assigned to more than one element
          if (obsrv%indxbnds_count > 1) then
            write (errmsg, '(a, a, a, a)') &
              trim(adjustl(this%depvartype))// &
              ' for observation', trim(adjustl(obsrv%Name)), &
              ' must be assigned to a feature with a unique boundname.'
            call store_error(errmsg)
          end if
        case ('LKT', 'SFT', 'MWT', 'UZT', 'LKE', 'SFE', 'MWE', 'UZE')
          call this%rp_obs_budterm(obsrv, &
                                   this%flowbudptr%budterm(this%idxbudgwf))
        case ('FLOW-JA-FACE')
          if (this%idxbudfjf > 0) then
            call this%rp_obs_flowjaface(obsrv, &
                                        this%flowbudptr%budterm(this%idxbudfjf))
          else
            write (errmsg, '(7a)') &
              'Observation ', trim(obsrv%Name), ' of type ', &
              trim(adjustl(obsrv%ObsTypeId)), ' in package ', &
              trim(this%packName), &
              ' cannot be processed because there are no flow connections.'
            call store_error(errmsg)
          end if
        case ('STORAGE')
          call this%rp_obs_byfeature(obsrv)
        case ('CONSTANT')
          call this%rp_obs_byfeature(obsrv)
        case ('FROM-MVR')
          call this%rp_obs_byfeature(obsrv)
        case default
          !
          ! -- check the child package for any specific obs
          found = .false.
          call this%pak_rp_obs(obsrv, found)
          !
          ! -- if none found then terminate with an error
          if (.not. found) then
            errmsg = 'Unrecognized observation type "'// &
                     trim(obsrv%ObsTypeId)//'" for '// &
                     trim(adjustl(this%text))//' package '// &
                     trim(this%packName)
            call store_error(errmsg, terminate=.TRUE.)
          end if
        end select

      end do
      !
      ! -- check for errors
      if (count_errors() > 0) then
        call store_error_unit(this%obs%inunitobs)
      end if
    end if
  end subroutine apt_rp_obs

  !> @brief Calculate observation values
  !!
  !! Routine calculates observations common to SFT/LKT/MWT/UZT
  !! (or SFE/LKE/MWE/UZE) for as many TspAptType observations that are common
  !! among the advanced transport packages
  !<
  subroutine apt_bd_obs(this)
    ! -- modules
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: igwfnode
    integer(I4B) :: j
    integer(I4B) :: jj
    integer(I4B) :: n
    integer(I4B) :: n1
    integer(I4B) :: n2
    real(DP) :: v
    type(ObserveType), pointer :: obsrv => null()
    logical :: found
    !
    ! -- Write simulated values for all Advanced Package observations
    if (this%obs%npakobs > 0) then
      call this%obs%obs_bd_clear()
      do i = 1, this%obs%npakobs
        obsrv => this%obs%pakobs(i)%obsrv
        do j = 1, obsrv%indxbnds_count
          v = DNODATA
          jj = obsrv%indxbnds(j)
          select case (obsrv%ObsTypeId)
          case ('CONCENTRATION', 'TEMPERATURE')
            if (this%iboundpak(jj) /= 0) then
              v = this%xnewpak(jj)
            end if
          case ('LKT', 'SFT', 'MWT', 'UZT', 'LKE', 'SFE', 'MWE', 'UZE')
            n = this%flowbudptr%budterm(this%idxbudgwf)%id1(jj)
            if (this%iboundpak(n) /= 0) then
              igwfnode = this%flowbudptr%budterm(this%idxbudgwf)%id2(jj)
              v = this%hcof(jj) * this%xnew(igwfnode) - this%rhs(jj)
              v = -v
            end if
          case ('FLOW-JA-FACE')
            n = this%flowbudptr%budterm(this%idxbudfjf)%id1(jj)
            if (this%iboundpak(n) /= 0) then
              call this%apt_fjf_term(jj, n1, n2, v)
            end if
          case ('STORAGE')
            if (this%iboundpak(jj) /= 0) then
              v = this%qsto(jj)
            end if
          case ('CONSTANT')
            if (this%iboundpak(jj) /= 0) then
              v = this%ccterm(jj)
            end if
          case ('FROM-MVR')
            if (this%iboundpak(jj) /= 0 .and. this%idxbudfmvr > 0) then
              call this%apt_fmvr_term(jj, n1, n2, v)
            end if
          case ('TO-MVR')
            if (this%idxbudtmvr > 0) then
              n = this%flowbudptr%budterm(this%idxbudtmvr)%id1(jj)
              if (this%iboundpak(n) /= 0) then
                call this%apt_tmvr_term(jj, n1, n2, v)
              end if
            end if
          case default
            found = .false.
            !
            ! -- check the child package for any specific obs
            call this%pak_bd_obs(obsrv%ObsTypeId, jj, v, found)
            !
            ! -- if none found then terminate with an error
            if (.not. found) then
              errmsg = 'Unrecognized observation type "'// &
                       trim(obsrv%ObsTypeId)//'" for '// &
                       trim(adjustl(this%text))//' package '// &
                       trim(this%packName)
              call store_error(errmsg, terminate=.TRUE.)
            end if
          end select
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      end do
      !
      ! -- write summary of error messages
      if (count_errors() > 0) then
        call store_error_unit(this%obs%inunitobs)
      end if
    end if
  end subroutine apt_bd_obs

  !> @brief Check if observation exists in an advanced package
  !<
  subroutine pak_bd_obs(this, obstypeid, jj, v, found)
    ! -- dummy
    class(TspAptType), intent(inout) :: this
    character(len=*), intent(in) :: obstypeid
    integer(I4B), intent(in) :: jj
    real(DP), intent(inout) :: v
    logical, intent(inout) :: found
    ! -- local
    !
    ! -- set found = .false. because obstypeid is not known
    found = .false.
  end subroutine pak_bd_obs

  !> @brief Process observation IDs for an advanced package
  !!
  !! Method to process observation ID strings for an APT package.
  !! This processor is only for observation types that support ID1
  !! and not ID2.
  !<
  subroutine apt_process_obsID(obsrv, dis, inunitobs, iout)
    ! -- dummy variables
    type(ObserveType), intent(inout) :: obsrv !< Observation object
    class(DisBaseType), intent(in) :: dis !< Discretization object
    integer(I4B), intent(in) :: inunitobs !< file unit number for the package observation file
    integer(I4B), intent(in) :: iout !< model listing file unit number
    ! -- local variables
    integer(I4B) :: nn1
    integer(I4B) :: icol
    integer(I4B) :: istart
    integer(I4B) :: istop
    character(len=LINELENGTH) :: string
    character(len=LENBOUNDNAME) :: bndname
    !
    ! -- initialize local variables
    string = obsrv%IDstring
    !
    ! -- Extract reach number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    boundary name--deal with it.
    icol = 1
    !
    ! -- get reach number or boundary name
    call extract_idnum_or_bndname(string, icol, istart, istop, nn1, bndname)
    if (nn1 == NAMEDBOUNDFLAG) then
      obsrv%FeatureName = bndname
    end if
    !
    ! -- store reach number (NodeNumber)
    obsrv%NodeNumber = nn1
    !
    ! -- store NodeNumber2 as 1 so that this can be used
    !    as the iconn value for SFT.  This works for SFT
    !    because there is only one reach per GWT connection.
    obsrv%NodeNumber2 = 1
  end subroutine apt_process_obsID

  !> @brief Process observation IDs for a package
  !!
  !! Method to process observation ID strings for an APT package. This
  !! processor is for the case where if ID1 is an integer then ID2 must be
  !! provided.
  !<
  subroutine apt_process_obsID12(obsrv, dis, inunitobs, iout)
    ! -- dummy variables
    type(ObserveType), intent(inout) :: obsrv !< Observation object
    class(DisBaseType), intent(in) :: dis !< Discretization object
    integer(I4B), intent(in) :: inunitobs !< file unit number for the package observation file
    integer(I4B), intent(in) :: iout !< model listing file unit number
    ! -- local variables
    integer(I4B) :: nn1
    integer(I4B) :: iconn
    integer(I4B) :: icol
    integer(I4B) :: istart
    integer(I4B) :: istop
    character(len=LINELENGTH) :: string
    character(len=LENBOUNDNAME) :: bndname
    !
    ! -- initialize local variables
    string = obsrv%IDstring
    !
    ! -- Extract reach number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    boundary name--deal with it.
    icol = 1
    !
    ! -- get reach number or boundary name
    call extract_idnum_or_bndname(string, icol, istart, istop, nn1, bndname)
    if (nn1 == NAMEDBOUNDFLAG) then
      obsrv%FeatureName = bndname
    else
      call extract_idnum_or_bndname(string, icol, istart, istop, iconn, bndname)
      if (len_trim(bndName) < 1 .and. iconn < 0) then
        write (errmsg, '(a,1x,a,a,1x,a,1x,a)') &
          'For observation type', trim(adjustl(obsrv%ObsTypeId)), &
          ', ID given as an integer and not as boundname,', &
          'but ID2 is missing.  Either change ID to valid', &
          'boundname or supply valid entry for ID2.'
        call store_error(errmsg)
      end if
      obsrv%NodeNumber2 = iconn
    end if
    !
    ! -- store reach number (NodeNumber)
    obsrv%NodeNumber = nn1
  end subroutine apt_process_obsID12

  !> @brief Setup a table object an advanced package
  !!
  !! Set up the table object that is used to write the apt concentration
  !! (or temperature) data. The terms listed here must correspond in the
  !! apt_ot method.
  !<
  subroutine apt_setup_tableobj(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH, LENBUDTXT
    ! -- dummy
    class(TspAptType) :: this
    ! -- local
    integer(I4B) :: nterms
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text_temp
    !
    ! -- setup well head table
    if (this%iprconc > 0) then
      !
      ! -- Determine the number of head table columns
      nterms = 2
      if (this%inamedbound == 1) nterms = nterms + 1
      !
      ! -- set up table title
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(adjustl(this%packName))// &
              ') '//trim(adjustl(this%depvartype))// &
              &' FOR EACH CONTROL VOLUME'
      !
      ! -- set up dv tableobj
      call table_cr(this%dvtab, this%packName, title)
      call this%dvtab%table_df(this%ncv, nterms, this%iout, &
                               transient=.TRUE.)
      !
      ! -- Go through and set up table budget term
      if (this%inamedbound == 1) then
        text_temp = 'NAME'
        call this%dvtab%initialize_column(text_temp, 20, alignment=TABLEFT)
      end if
      !
      ! -- feature number
      text_temp = 'NUMBER'
      call this%dvtab%initialize_column(text_temp, 10, alignment=TABCENTER)
      !
      ! -- feature conc
      text_temp = this%depvartype(1:4)
      call this%dvtab%initialize_column(text_temp, 12, alignment=TABCENTER)
    end if
  end subroutine apt_setup_tableobj

end module TspAptModule
