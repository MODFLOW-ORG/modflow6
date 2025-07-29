!GWF Water Mover Module
!This module contains a derived type, called GwfMvrType, that
!is attached to the GWF model.  The water mover can be used to move water
!between packages.  The mover requires that mover-aware packages have access
!to four arrays: qtformvr, qformvr, qtomvr, and qfrommvr.  These arrays are
!stored and managed by a separate PackageMoverType object.  qformvr is a
!vector of volumetric flow rates available for the mover.  The package
!must fill the vector (dimensioned by number of reaches) with the available
!water.  qtomvr is a vector containing how much water was actually moved
!by the mover.  The package should use this value in the budgeting part
!to track how much water was actually provided to the mover.  Lastly,
!the qfrommvr is a vector that contains volumetric rates for how much
!water was provided by the mover as a source of water to the package.
!
!The mover is designed so that a reach can provide water to more than one
!receiving reaches.  The available water will be consumed in order of
!the movers listed in the package.  The mover is also designed so that
!a receiver can receive water from more than one provider.
!
!  1.  The mover is instantiated as a model member:
!
!      type(GwfMvrType),               pointer :: mvr     => null()
!
!      Mover aware packages have access to the following vectors of mover
!      information, which are stored in the PackageMoverType object:
!
!      qtformvr(nproviders) -- total available unconsumed water for mover
!      qformvr(nproviders) -- currently available consumed water (changes during fc)
!      qtomvr(nproviders) -- actual amount of water sent to mover
!      qfrommvr(nreceivers) -- actual amount of water received from mover
!
!      integer(I4B), pointer                        :: imover        => null()
!      real(DP), dimension(:), pointer, contiguous  :: qtformvr      => null()
!      real(DP), dimension(:), pointer, contiguous  :: qformvr       => null()
!      real(DP), dimension(:), pointer, contiguous  :: qtomvr        => null()
!      real(DP), dimension(:), pointer, contiguous  :: qfrommvr      => null()
!
!      Note qtformvr is filled as a positive number to indicate that it is
!      water available to be moved.  If qtformvr is negative, then
!      no water will be moved for that reach.  qformvr is also the available
!      water, but this value decreases as the mover object consumes water from
!      it.
!
!  2.  In gwf_cr create the mover package by calling the CR subroutine:
!
!      call mvr_cr(this%mvr, this%name, this%inmvr, this%iout)
!
!  3.  In gwf_ar call the AR method for the mover:
!
!      if(this%inmvr > 0) call this%mvr%mvr_ar()
!
!      Mover aware packages allocate the four vectors.  The first three
!      (qtformvr, qformvr, qtomvr) are allocated to the number of providers
!      and the last one (qfrommvr) is allocated to the number of receivers.
!
!  4.  In gwf_rp call the RP method for the mover.  This reads the
!      movers active for the current period.
!
!      if(this%inmvr > 0) call this%mvr%mvr_rp()
!
!  5.  In gwf_ad call the AD method for the mover.  This saves qtomvr from the
!      the last time step.
!
!      if(this%inmvr > 0) call this%mvr%mvr_ad()
!
!      Mover aware packages then set:
!        qtomvr(:) = 0.
!        qformvr(:) = 0.
!
!  6.  In gwf_cf call the CF routine. Mover aware packages set:
!        qtformvr(:) = qformvr(:)
!        qfrommvr(:) = 0.
!        qtomvr(:) = 0.
!
!  7.  The FC method for the mover is called.  This method calculates the
!      amount of water to move based on the amount of water available from the
!      previous iteration.  This call updates the values in the qtomvr and
!      qfrommvr vectors inside the packages.  This is done by the mover package
!      using pointers to the appropriate reach locations in qtomvr and qfrommvr.
!
!      if(this%inmvr > 0) call this%mvr%mvr_fc()  ! called from gwf%gwf_fc()
!
!      a. Mover aware packages first set qformvr(:) = 0.
!      b. Mover aware packages that are receivers (MAW, SFR, LAK, UZF) add
!         qfrommvr terms to their individual control volume equations as a
!         source of water.
!      c. Mover aware packages calculate qformvr as amount of water available
!         to be moved (these qformvr terms are used in the next iteration
!         by this%mvr%mvr_fc() to calculate how much water is actually moved)
!
!  8.  The BD method for the mover is called.  This method writes the moved
!      water rates if requested.
!
!      if(this%inmvr > 0) call this%mvr%mvr_bd()
!
!      Mover aware packages account for qtomvr and qfrommvr terms in their
!      individual budget routines.
!
!  9.  The OT method for the mover is called.  This method outputs a mover
!      budget table.
!
!      if(this%inmvr > 0) call this%mvr%mvr_ot()
!
module GwfMvrModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: LENMEMPATH, LENPACKAGENAME, LENMODELNAME, &
                             LENBUDTXT, LENAUXNAME, LENPAKLOC, &
                             DZERO, DNODATA, MAXCHARLEN, TABCENTER, &
                             LINELENGTH
  use SimModule, only: store_error, store_error_filename, count_errors
  use MvrModule, only: MvrType
  use BudgetModule, only: BudgetType, budget_cr
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr
  use NumericalPackageModule, only: NumericalPackageType
  use PackageMoverModule, only: PackageMoverType, set_packagemover_pointer
  use BaseDisModule, only: DisBaseType
  use InputOutputModule, only: urword
  use TableModule, only: TableType, table_cr

  implicit none
  private
  public :: GwfMvrType, mvr_cr

  type, extends(NumericalPackageType) :: GwfMvrType
    logical(LGP), pointer :: reset_mapped_id ! flag to indicate mapped ids must be reset; true when movers change
    integer(I4B), pointer :: ibudgetout => null() !< binary budget output file
    integer(I4B), pointer :: ibudcsv => null() !< unit number for csv budget output file
    integer(I4B), pointer :: maxmvr => null() !< max number of movers to be specified
    integer(I4B), pointer :: maxpackages => null() !< max number of packages to be specified
    integer(I4B), pointer :: maxcomb => null() !< max number of combination of packages
    integer(I4B), pointer :: nmvr => null() !< number of movers for current stress period
    integer(I4B), pointer :: iexgmvr => null() !< indicate mover is for an exchange (not for a single model)
    integer(I4B), pointer :: imodelnames => null() !< indicate package input file has model names in it
    integer(I4B), dimension(:), pointer, contiguous :: ientries => null() !< number of entries for each combination
    character(len=LENMEMPATH), &
      dimension(:), pointer, contiguous :: pckMemPaths !< memory paths of all packages used in this mover
    character(len=LENPACKAGENAME), &
      dimension(:), pointer, contiguous :: paknames => null() !< array of package names
    type(MvrType), dimension(:), pointer, contiguous :: mvr => null() !< array of movers
    type(BudgetType), pointer :: budget => null() !< mover budget object (used to write table)
    type(BudgetObjectType), pointer :: budobj => null() !< new budget container (used to write binary file)
    type(PackageMoverType), &
      dimension(:), pointer, contiguous :: pakmovers => null() !< pointer to package mover objects
    !
    ! -- table objects
    type(TableType), pointer :: outputtab => null()
    logical(LGP) :: suppress_fileout = .false. !< flag to disable output file (budget, budget csv)

  contains
    procedure :: mvr_init
    procedure :: mvr_ar
    procedure :: mvr_rp
    procedure :: mvr_ad
    procedure :: mvr_fc
    procedure :: mvr_cc
    procedure :: mvr_bd
    procedure :: mvr_bdsav
    procedure :: mvr_ot_saveflow
    procedure :: mvr_ot_printflow
    procedure :: mvr_ot_bdsummary
    procedure :: mvr_da
    procedure :: source_options
    procedure :: check_options
    procedure :: log_options
    procedure :: source_dimensions
    procedure :: source_packages
    procedure :: check_packages
    procedure :: assign_packagemovers
    procedure :: initialize_movers
    procedure :: fill_budobj
    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure, private :: mvr_setup_budobj
    procedure, private :: mvr_setup_outputtab
    procedure, private :: mvr_print_outputtab
    procedure, private :: set_mapped_id

  end type GwfMvrType

contains

  !> @brief Create a new mvr object
  !<
  subroutine mvr_cr(mvrobj, name_parent, input_mempath, inunit, iout, dis, &
                    iexgmvr)
    ! -- dummy
    type(GwfMvrType), pointer :: mvrobj
    character(len=*), intent(in) :: name_parent
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), optional :: iexgmvr
    !
    ! -- Create the object
    allocate (mvrobj)
    !
    ! -- Init
    call mvrobj%mvr_init(name_parent, input_mempath, inunit, iout, dis, iexgmvr)
  end subroutine mvr_cr

  subroutine mvr_init(this, name_parent, input_mempath, inunit, iout, dis, &
                      iexgmvr)
    class(GwfMvrType) :: this
    character(len=*), intent(in) :: name_parent
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    class(DisBaseType), pointer, intent(in) :: dis
    integer(I4B), optional :: iexgmvr
    !
    ! -- create name and memory paths. name_parent will either be model name or the
    !    exchange name.
    call this%set_names(1, name_parent, 'MVR', 'MVR', input_mempath)
    !
    ! -- Allocate scalars
    call this%allocate_scalars()
    !
    ! -- Set pointer to dis
    this%dis => dis
    !
    ! -- Set variables
    this%inunit = inunit
    this%iout = iout
    !
    ! -- Set iexgmvr
    if (present(iexgmvr)) this%iexgmvr = iexgmvr
    !
    ! -- Create the budget object
    if (inunit > 0) then
      call budget_cr(this%budget, this%memoryPath)
    end if
    !
    ! -- instantiate the budget object
    call budgetobject_cr(this%budobj, 'WATER MOVER')
  end subroutine mvr_init

  !> @brief Allocate and read water mover information
  !<
  subroutine mvr_ar(this)
    ! -- dummy
    class(GwfMvrType) :: this
    !
    ! -- Print a message identifying the water mover package.
    write (this%iout, 1) this%input_mempath
1   format(1x, /1x, 'MVR -- WATER MOVER PACKAGE, VERSION 8, 1/29/2016', &
           ' INPUT READ FROM MEMPATH: ', A)
    !
    ! -- Read and check options
    call this%source_options()
    call this%check_options()
    !
    ! -- Read options
    call this%source_dimensions()
    !
    ! -- Allocate arrays
    call this%allocate_arrays()
    !
    ! -- Read and check package names
    call this%source_packages()
    call this%check_packages()
    !
    ! -- Define the budget object to be the size of package names
    call this%budget%budget_df(this%maxpackages, 'WATER MOVER')
    call this%budget%set_ibudcsv(this%ibudcsv)
    !
    ! -- setup the budget object
    call this%mvr_setup_budobj()
  end subroutine mvr_ar

  !> @brief Read and Prepare
  !!
  !! Read itmp and read new boundaries if itmp > 0
  !<
  subroutine mvr_rp(this)
    ! -- modules
    use TdisModule, only: kper
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use ArrayHandlersModule, only: ifind
    ! -- dummy
    class(GwfMvrType), intent(inout) :: this
    ! -- local
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: mnames1, mnames2, pnames1, pnames2, mvrtypes
    integer(I4B), dimension(:), pointer, contiguous :: id1, id2
    real(DP), dimension(:), pointer, contiguous :: values
    character(len=LINELENGTH) :: mname1, mname2, pname1, pname2, mvrtype
    integer(I4B), pointer :: iper
    integer(I4B) :: n, ii, jj, imvrtype, ipos
    logical(LGP) :: found
    character(len=LINELENGTH) :: errmsg
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Error.  Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    character(len=*), parameter :: fmtnbd = &
      "(1X,/1X,'THE NUMBER OF ACTIVE ',A,'S (',I6, &
       &') IS GREATER THAN MAXIMUM(',I6,')')"

    ! -- Assign a pointer to the package mover object.  The pointer assignment
    !    will happen only the first time
    call this%assign_packagemovers()

    call mem_setptr(iper, 'IPER', this%input_mempath)
    if (iper /= kper) then
      write (this%iout, fmtlsp) 'MVR'
      return
    end if

    ! force mapped id reset
    this%reset_mapped_id = .true.

    ! reset and update nmvr
    this%nmvr = 0
    call mem_set_value(this%nmvr, 'NBOUND', this%input_mempath, &
                       found)

    ! set ponters to input context
    if (this%imodelnames /= 0) then
      call mem_setptr(mnames1, 'MNAME1', this%input_mempath)
      call mem_setptr(mnames2, 'MNAME2', this%input_mempath)
    end if
    call mem_setptr(pnames1, 'PNAME1', this%input_mempath)
    call mem_setptr(pnames2, 'PNAME2', this%input_mempath)
    call mem_setptr(id1, 'ID1', this%input_mempath)
    call mem_setptr(id2, 'ID2', this%input_mempath)
    call mem_setptr(mvrtypes, 'MVRTYPE', this%input_mempath)
    call mem_setptr(values, 'VALUE', this%input_mempath)

    write (this%iout, '(/,1x,a,i0)') 'READING WATER MOVERS FOR PERIOD ', kper
    !
    ! process input
    do n = 1, this%nmvr
      if (this%imodelnames /= 0) then
        mname1 = mnames1(n)
        mname2 = mnames2(n)
      else
        mname1 = this%name_model
        mname2 = this%name_model
      end if
      pname1 = pnames1(n)
      pname2 = pnames2(n)
      mvrtype = mvrtypes(n)

      select case (mvrtype)
      case ('FACTOR')
        imvrtype = 1
      case ('EXCESS')
        imvrtype = 2
      case ('THRESHOLD')
        imvrtype = 3
      case ('UPTO')
        imvrtype = 4
      case default
        imvrtype = 0
      end select

      call this%mvr(n)%set_values(mname1, pname1, id1(n), &
                                  mname2, pname2, id2(n), &
                                  imvrtype, values(n))
    end do

    ! Process the input data into the individual mover objects
    call this%initialize_movers(this%nmvr)

    ! assign the pointers
    do n = 1, this%nmvr
      call this%mvr(n)%prepare(this%pckMemPaths, &
                               this%pakmovers, &
                               this%input_fname)
      if (this%iprpak == 1) call this%mvr(n)%echo(this%iout)
    end do

    ! Check to make sure all providers and receivers are properly stored
    do n = 1, this%nmvr
      ipos = ifind(this%pckMemPaths, this%mvr(n)%mem_path_src)
      if (ipos < 1) then
        write (errmsg, '(a,a,a)') 'Provider ', &
          trim(this%mvr(n)%mem_path_src), ' not listed in packages block.'
        call store_error(errmsg)
      end if
      ipos = ifind(this%pckMemPaths, this%mvr(n)%mem_path_tgt)
      if (ipos < 1) then
        write (errmsg, '(a,a,a)') 'Receiver ', &
          trim(this%mvr(n)%mem_path_tgt), ' not listed in packages block.'
        call store_error(errmsg)
      end if
    end do
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    ! reset ientries
    do n = 1, this%maxcomb
      this%ientries(n) = 0
    end do

    ! set ientries
    do n = 1, this%nmvr
      ii = ifind(this%pckMemPaths, this%mvr(n)%mem_path_src)
      jj = ifind(this%pckMemPaths, this%mvr(n)%mem_path_tgt)
      ipos = (ii - 1) * this%maxpackages + jj
      this%ientries(ipos) = this%ientries(ipos) + 1
      write (this%iout, '(3x,a,i0)') 'IENTRY ipos=', ipos
    end do

    write (this%iout, '(1x,a,1x,i6)') 'END OF INPUT FOR PERIOD', kper
    write (this%iout, '(1x, i0, a, i0,/)') this%nmvr, &
      ' MOVERS READ FOR PERIOD ', kper
  end subroutine mvr_rp

  subroutine initialize_movers(this, nr_active_movers)
    class(GwfMvrType) :: this
    integer(I4B) :: nr_active_movers
    ! local
  end subroutine initialize_movers

  subroutine mvr_ad(this)
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- locals
    integer(I4B) :: i
    !
    do i = 1, this%nmvr
      call this%mvr(i)%advance()
    end do
  end subroutine mvr_ad

  !> @brief Calculate qfrommvr as a function of qtomvr
  !<
  subroutine mvr_fc(this)
    class(GwfMvrType) :: this
    ! local
    integer(I4B) :: i

    do i = 1, this%nmvr
      call this%mvr(i)%update_provider()
      call this%mvr(i)%update_receiver()
    end do

  end subroutine mvr_fc

  !> @brief Extra convergence check for mover
  !<
  subroutine mvr_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- dummy
    class(GwfMvrType) :: this
    integer(I4B), intent(in) :: innertot
    integer(I4B), intent(in) :: kiter
    integer(I4B), intent(in) :: iend
    integer(I4B), intent(in) :: icnvgmod
    character(len=LENPAKLOC), intent(inout) :: cpak
    integer(I4B), intent(inout) :: ipak
    real(DP), intent(inout) :: dpak
    ! -- formats
    character(len=*), parameter :: fmtmvrcnvg = &
      "(/,1x,'MOVER PACKAGE REQUIRES AT LEAST TWO OUTER ITERATIONS. CONVERGE &
      &FLAG HAS BEEN RESET TO FALSE.')"
    !
    ! -- If there are active movers, then at least 2 outers required
    if (this%nmvr > 0) then
      if (icnvgmod == 1 .and. kiter == 1) then
        dpak = DNODATA
        cpak = trim(this%packName)
        write (this%iout, fmtmvrcnvg)
      end if
    end if
  end subroutine mvr_cc

  !> @brief Fill the mover budget object
  !<
  subroutine mvr_bd(this)
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- locals
    ! -- formats
    !
    ! -- set the feature maps; for performance reasons,
    !    this should only be called for the first time
    !    step of a stress period in which a new set of
    !    movers was provided in a period block.
    if (this%reset_mapped_id) then
      call this%set_mapped_id()
      this%reset_mapped_id = .false.
    end if
    !
    ! -- fill the budget object
    call this%fill_budobj()
  end subroutine mvr_bd

  !> @brief Write mover terms
  !<
  subroutine mvr_bdsav(this, icbcfl, ibudfl, isuppress_output)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, pertim, totim
    ! -- dummy
    class(GwfMvrType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    integer(I4B), intent(in) :: isuppress_output
    ! -- locals
    integer(I4B) :: ibinun
    ! -- formats
    character(len=*), parameter :: fmttkk = &
                                   "(1X,/1X,A,'   PERIOD ',I0,'   STEP ',I0)"
    !
    ! -- Print the mover flow table
    if (ibudfl /= 0 .and. this%iprflow /= 0 .and. isuppress_output == 0) then
      call this%mvr_print_outputtab()
    end if
    !
    ! -- Save the mover flows from the budobj to a mover binary file
    ibinun = 0
    if (this%ibudgetout /= 0) then
      ibinun = this%ibudgetout
    end if
    if (icbcfl == 0) ibinun = 0
    if (isuppress_output /= 0) ibinun = 0
    if (ibinun > 0) then
      call this%budobj%save_flows(this%dis, ibinun, kstp, kper, delt, &
                                  pertim, totim, this%iout)
    end if
  end subroutine mvr_bdsav

  !> @brief Write mover terms
  !<
  subroutine mvr_ot_saveflow(this, icbcfl, ibudfl)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, pertim, totim
    ! -- dummy
    class(GwfMvrType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    ! -- locals
    integer(I4B) :: ibinun
    !
    ! -- Save the mover flows from the budobj to a mover binary file
    ibinun = 0
    if (this%ibudgetout /= 0) then
      ibinun = this%ibudgetout
    end if
    if (icbcfl == 0) ibinun = 0
    if (ibinun > 0) then
      call this%budobj%save_flows(this%dis, ibinun, kstp, kper, delt, &
                                  pertim, totim, this%iout)
    end if
  end subroutine mvr_ot_saveflow

  !> @brief Print mover flow table
  !<
  subroutine mvr_ot_printflow(this, icbcfl, ibudfl)
    ! -- dummy
    class(GwfMvrType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: ibudfl
    !
    ! -- Print the mover flow table
    if (ibudfl /= 0 .and. this%iprflow /= 0) then
      call this%mvr_print_outputtab()
    end if
  end subroutine mvr_ot_printflow

  !> @brief Write mover budget to listing file
  !<
  subroutine mvr_ot_bdsummary(this, ibudfl)
    ! -- modules
    use TdisModule, only: kstp, kper, delt, totim
    ! -- dummy
    class(GwfMvrType) :: this
    integer(I4B), intent(in) :: ibudfl
    ! -- locals
    character(len=LENMEMPATH) :: pckMemPath
    integer(I4B) :: i, j
    real(DP), allocatable, dimension(:) :: ratin, ratout
    !
    ! -- Allocate and initialize ratin/ratout
    allocate (ratin(this%maxpackages), ratout(this%maxpackages))
    do j = 1, this%maxpackages
      ratin(j) = DZERO
      ratout(j) = DZERO
    end do
    !
    ! -- Accumulate the rates
    do i = 1, this%nmvr
      do j = 1, this%maxpackages
        if (this%pckMemPaths(j) == this%mvr(i)%mem_path_src) then
          ratin(j) = ratin(j) + this%mvr(i)%qpactual
        end if
        if (this%pckMemPaths(j) == this%mvr(i)%mem_path_tgt) then
          ratout(j) = ratout(j) + this%mvr(i)%qpactual
        end if
      end do
    end do
    !
    ! -- Send rates to budget object
    call this%budget%reset()
    do j = 1, this%maxpackages
      if ((this%iexgmvr) == 1) then
        pckMemPath = this%pckMemPaths(j)
      else
        pckMemPath = this%paknames(j)
      end if
      call this%budget%addentry(ratin(j), ratout(j), delt, pckMemPath)
    end do
    !
    ! -- Write the budget
    call this%budget%finalize_step(delt)
    if (ibudfl /= 0) then
      call this%budget%budget_ot(kstp, kper, this%iout)
    end if
    !
    ! -- Write budget csv
    call this%budget%writecsv(totim)
    !
    ! -- Deallocate
    deallocate (ratin, ratout)
    !
    ! -- Output mvr budget
    !    Not using budobj write_table here because it would result
    !    in a table that has one entry.  A custom table looks
    !    better here with a row for each package.
    !call this%budobj%write_budtable(kstp, kper, this%iout)
  end subroutine mvr_ot_bdsummary

  !> @brief Deallocate
  !<
  subroutine mvr_da(this)
    ! -- modules
    use ConstantsModule, only: DONE
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfMvrType) :: this
    !
    ! -- Arrays
    if (this%inunit > 0) then
      call mem_deallocate(this%ientries)
      deallocate (this%mvr)
      deallocate (this%pckMemPaths)
      deallocate (this%paknames)
      deallocate (this%pakmovers)
      !
      ! -- budget object
      call this%budget%budget_da()
      deallocate (this%budget)
      !
      ! -- budobj
      call this%budobj%budgetobject_da()
      deallocate (this%budobj)
      nullify (this%budobj)
      !
      ! -- output table object
      if (associated(this%outputtab)) then
        call this%outputtab%table_da()
        deallocate (this%outputtab)
        nullify (this%outputtab)
      end if
    end if
    !
    ! -- Scalars
    call mem_deallocate(this%reset_mapped_id)
    call mem_deallocate(this%ibudgetout)
    call mem_deallocate(this%ibudcsv)
    call mem_deallocate(this%maxmvr)
    call mem_deallocate(this%maxpackages)
    call mem_deallocate(this%maxcomb)
    call mem_deallocate(this%nmvr)
    call mem_deallocate(this%iexgmvr)
    call mem_deallocate(this%imodelnames)
    !
    ! -- deallocate scalars in NumericalPackageType
    call this%NumericalPackageType%da()
  end subroutine mvr_da

  subroutine source_options(this)
    ! -- modules
    use OpenSpecModule, only: access, form
    use InputOutputModule, only: assign_iounit, openfile
    use MemoryManagerExtModule, only: mem_set_value
    use GwfMvrInputModule, only: GwfMvrParamFoundType
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    type(GwfMvrParamFoundType) :: found
    character(len=LINELENGTH) :: budgetfile, budgetcsvfile

    !
    ! -- source package input
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, &
                       found%print_input)
    call mem_set_value(this%iprflow, 'PRINT_FLOWS', this%input_mempath, &
                       found%print_flows)
    call mem_set_value(this%imodelnames, 'MODELNAMES', this%input_mempath, &
                       found%modelnames)
    call mem_set_value(budgetfile, 'BUDGETFILE', this%input_mempath, &
                       found%budgetfile)
    call mem_set_value(budgetcsvfile, 'BUDGETCSVFILE', this%input_mempath, &
                       found%budgetcsvfile)

    ! modelnames is validated in check_options()

    ! fileout options
    if (found%budgetfile) then
      if (.not. this%suppress_fileout) then
        call assign_iounit(this%ibudgetout, this%inunit, "BUDGET fileout")
        call openfile(this%ibudgetout, this%iout, budgetfile, 'DATA(BINARY)', &
                      form, access, 'REPLACE')
      end if
    end if
    if (found%budgetcsvfile) then
      if (.not. this%suppress_fileout) then
        call assign_iounit(this%ibudcsv, this%inunit, "BUDGETCSV fileout")
        call openfile(this%ibudcsv, this%iout, budgetcsvfile, 'CSV', &
                      filstat_opt='REPLACE')
      end if
    end if

    call this%log_options(found, budgetfile, budgetcsvfile)
  end subroutine source_options

  !> @brief Log MVR options
  !<
  subroutine log_options(this, found, budgetfile, budgetcsvfile)
    ! -- modules
    use GwfMvrInputModule, only: GwfMvrParamFoundType
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    type(GwfMvrParamFoundType), intent(in) :: found
    character(len=*), intent(in) :: budgetfile
    character(len=*), intent(in) :: budgetcsvfile
    ! -- formats
    character(len=*), parameter :: fmtmvrbin = &
      "(4x, 'MVR ', 1x, a, 1x, ' WILL BE SAVED TO FILE: ', a, /4x, 'OPENED ON &
      &UNIT: ', I0)"

    write (this%iout, '(1x,a)') 'PROCESSING MVR OPTIONS'

    if (found%budgetfile) then
      if (.not. this%suppress_fileout) then
        write (this%iout, fmtmvrbin) 'BUDGET', trim(adjustl(budgetfile)), &
          this%ibudgetout
      end if
    end if

    if (found%budgetcsvfile) then
      if (.not. this%suppress_fileout) then
        write (this%iout, fmtmvrbin) 'BUDGET CSV', trim(adjustl(budgetcsvfile)), &
          this%ibudcsv
      end if
    end if

    if (found%print_input) then
      write (this%iout, '(4x,a)') 'WATER MOVER INPUT '// &
        'WILL BE PRINTED TO LIST FILE.'
    end if

    if (found%print_flows) then
      write (this%iout, '(4x,a)') 'LISTS OF WATER MOVER FLOWS '// &
        'WILL BE PRINTED TO LIST FILE.'
    end if

    if (found%modelnames) then
      if (this%iexgmvr /= 0) then
        write (this%iout, '(4x,a)') 'ALL PACKAGE NAMES ARE PRECEDED '// &
          'BY THE NAME OF THE MODEL CONTAINING THE PACKAGE.'
      end if
    end if

    write (this%iout, '(1x,a)') 'END OF MVR OPTIONS'
  end subroutine log_options

  !> @brief Check MODELNAMES option set correctly
  !<
  subroutine check_options(this)
    ! -- modules
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg
    !
    ! -- Check if not exchange mover but model names are specified
    if (this%iexgmvr == 0 .and. this%imodelnames == 1) then
      write (errmsg, '(a,a)') &
        'MODELNAMES cannot be specified unless the '// &
        'mover package is for an exchange.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- Check if exchange mover but model names not specified
    if (this%iexgmvr /= 0 .and. this%imodelnames == 0) then
      write (errmsg, '(a,a)') &
        'MODELNAMES option must be specified because '// &
        'mover package is for an exchange.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
  end subroutine check_options

  subroutine source_dimensions(this)
    ! -- modules
    use SimVariablesModule, only: errmsg
    use MemoryManagerExtModule, only: mem_set_value
    use GwfMvrInputModule, only: GwfMvrParamFoundType
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    type(GwfMvrParamFoundType) :: found
    integer(I4B) :: i, j

    call mem_set_value(this%maxmvr, 'MAXBOUND', this%input_mempath, &
                       found%maxbound)
    call mem_set_value(this%maxpackages, 'NPACKAGES', this%input_mempath, &
                       found%npackages)

    write (this%iout, '(/1x,a)') 'PROCESSING MVR DIMENSIONS'
    write (this%iout, '(4x,a,i0)') 'MAXMVR = ', this%maxmvr
    write (this%iout, '(4x,a,i0)') 'MAXPACKAGES = ', this%maxpackages
    write (this%iout, '(1x,a)') 'END OF MVR DIMENSIONS'

    ! -- verify dimensions were set
    if (this%maxmvr < 0) then
      write (errmsg, '(a)') &
        'MAXMVR was not specified or was specified incorrectly.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if
    if (this%maxpackages < 0) then
      write (errmsg, '(a)') &
        'MAXPACKAGES was not specified or was specified incorrectly.'
      call store_error(errmsg)
      call store_error_filename(this%input_fname)
    end if

    ! -- calculate maximum number of combinations
    this%maxcomb = 0
    do i = 1, this%maxpackages
      do j = 1, this%maxpackages
        this%maxcomb = this%maxcomb + 1
      end do
    end do
  end subroutine source_dimensions

  subroutine source_packages(this)
    ! -- modules
    use MemoryHelperModule, only: create_mem_path
    use MemoryManagerModule, only: mem_setptr
    use CharacterStringModule, only: CharacterStringType
    ! -- dummy
    class(GwfMvrType), intent(inout) :: this
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: mnames, pnames
    character(len=LINELENGTH) :: mname, pname, errmsg
    integer(I4B) :: n

    ! -- set input context pointers
    if (this%imodelnames /= 0) then
      call mem_setptr(mnames, 'MNAME', this%input_mempath)
    end if
    call mem_setptr(pnames, 'PNAME', this%input_mempath)

    if (size(pnames) /= this%maxpackages) then
      write (errmsg, '(a,i0,a)') 'ERROR: PACKAGES block number of rows (', &
        size(pnames), ') is not equivalent to dimension MAXPACKAGES.'
      call store_error_filename(this%input_fname)
    end if

    write (this%iout, '(/1x,a)') 'PROCESSING MVR PACKAGES'

    do n = 1, size(pnames)
      if (this%imodelnames /= 0) mname = mnames(n)
      pname = pnames(n)

      if (this%imodelnames /= 0) then
        this%pckMemPaths(n) = create_mem_path(mname, pname)
      else
        this%pckMemPaths(n) = create_mem_path(this%name_model, pname)
      end if

      this%paknames(n) = trim(pname)

      write (this%iout, '(3x,a,a)') 'INCLUDING PACKAGE: ', &
        trim(this%pckMemPaths(n))
    end do

    write (this%iout, '(1x,a)') 'END OF MVR PACKAGES'
  end subroutine source_packages

  !> @brief Check to make sure packages have mover activated
  !<
  subroutine check_packages(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(GwfMvrType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: errmsg
    integer(I4B) :: i
    integer(I4B), pointer :: imover_ptr
    !
    ! -- Check to make sure mover is activated for each package
    do i = 1, size(this%pckMemPaths)
      imover_ptr => null()
      call mem_setptr(imover_ptr, 'IMOVER', trim(this%pckMemPaths(i)))
      if (imover_ptr == 0) then
        write (errmsg, '(a, a, a)') &
          'ERROR.  MODEL AND PACKAGE "', &
          trim(this%pckMemPaths(i)), &
          '" DOES NOT HAVE MOVER SPECIFIED IN OPTIONS BLOCK.'
        call store_error(errmsg)
      end if
    end do
    !
    ! -- Terminate if errors detected.
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine check_packages

  !> @brief Assign pointer to each package's packagemover object
  !<
  subroutine assign_packagemovers(this)
    ! -- modules
    use PackageMoverModule, only: set_packagemover_pointer
    ! -- dummy
    class(GwfMvrType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i
    !
    ! -- Assign the package mover pointer if it hasn't been assigned yet
    do i = 1, size(this%pckMemPaths)
      if (this%pakmovers(i)%memoryPath == '') then
        call set_packagemover_pointer(this%pakmovers(i), &
                                      trim(this%pckMemPaths(i)))
      end if
    end do
  end subroutine assign_packagemovers

  !> @brief Allocate package scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use ConstantsModule, only: DONE
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfMvrType) :: this
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%reset_mapped_id, 'RESET_MAPPED_ID', this%memoryPath)
    call mem_allocate(this%ibudgetout, 'IBUDGETOUT', this%memoryPath)
    call mem_allocate(this%ibudcsv, 'IBUDCSV', this%memoryPath)
    call mem_allocate(this%maxmvr, 'MAXMVR', this%memoryPath)
    call mem_allocate(this%maxpackages, 'MAXPACKAGES', this%memoryPath)
    call mem_allocate(this%maxcomb, 'MAXCOMB', this%memoryPath)
    call mem_allocate(this%nmvr, 'NMVR', this%memoryPath)
    call mem_allocate(this%iexgmvr, 'IEXGMVR', this%memoryPath)
    call mem_allocate(this%imodelnames, 'IMODELNAMES', this%memoryPath)
    !
    ! -- Initialize
    this%reset_mapped_id = .false.
    this%ibudgetout = 0
    this%ibudcsv = 0
    this%maxmvr = -1
    this%maxpackages = -1
    this%maxcomb = 0
    this%nmvr = 0
    this%iexgmvr = 0
    this%imodelnames = 0
  end subroutine allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    use ConstantsModule, only: DZERO
    use PackageMoverModule, only: nulllify_packagemover_pointer
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    integer(I4B) :: i
    !
    ! -- Allocate
    allocate (this%mvr(this%maxmvr))
    allocate (this%pckMemPaths(this%maxpackages))
    allocate (this%paknames(this%maxpackages))
    allocate (this%pakmovers(this%maxpackages))
    !
    ! -- nullify the pakmovers
    do i = 1, this%maxpackages
      call nulllify_packagemover_pointer(this%pakmovers(i))
    end do
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%ientries, this%maxcomb, 'IENTRIES', this%memoryPath)
    !
    ! -- initialize ientries
    do i = 1, this%maxcomb
      this%ientries(i) = 0
    end do
    !
    ! -- setup the output table
    call this%mvr_setup_outputtab()
  end subroutine allocate_arrays

  !> @brief Set up the budget object that stores all the mvr flows
  !<
  subroutine mvr_setup_budobj(this)
    ! -- modules
    use ConstantsModule, only: LENBUDTXT
    use MemoryHelperModule, only: split_mem_path
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    integer(I4B) :: nbudterm
    integer(I4B) :: ncv
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: naux
    character(len=LENMODELNAME) :: modelname1, modelname2
    character(len=LENPACKAGENAME) :: packagename1, packagename2
    integer(I4B) :: maxlist
    integer(I4B) :: idx
    character(len=LENBUDTXT) :: text
    !
    ! -- Determine the number of mover budget terms. These are fixed for
    !    the simulation and cannot change.  A separate term is required
    !    for each possible provider/receiver combination.
    nbudterm = 0
    do i = 1, this%maxpackages
      do j = 1, this%maxpackages
        nbudterm = nbudterm + 1
      end do
    end do
    !
    ! -- Number of control volumes is set to be 0, because there aren't
    !    any for the mover
    ncv = 0
    !
    ! -- set up budobj
    call this%budobj%budgetobject_df(ncv, nbudterm, 0, 0)
    idx = 0
    !
    ! -- Go through and set up each budget term
    text = '      MOVER-FLOW'
    maxlist = this%maxmvr
    naux = 0
    do i = 1, this%maxpackages

      call split_mem_path(this%pckMemPaths(i), modelname1, packagename1)

      do j = 1, this%maxpackages

        idx = idx + 1
        call split_mem_path(this%pckMemPaths(j), modelname2, packagename2)
        call this%budobj%budterm(idx)%initialize(text, &
                                                 modelname1, &
                                                 packagename1, &
                                                 modelname2, &
                                                 packagename2, &
                                                 maxlist, .false., .false., &
                                                 naux)
      end do
    end do
  end subroutine mvr_setup_budobj

  !> @brief Fill budget object
  !<
  subroutine fill_budobj(this)
    ! -- modules
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- local
    integer(I4B) :: idx
    integer(I4B) :: i
    integer(I4B) :: j
    integer(I4B) :: n, n1, n2
    integer(I4B) :: ipos
    integer(I4B) :: ival
    integer(I4B) :: nitems
    integer(I4B) :: lloc
    integer(I4B) :: istart
    integer(I4B) :: istop
    real(DP) :: rval
    character(len=LENMODELNAME) :: modelname1, modelname2
    character(len=LENPACKAGENAME) :: packagename1, packagename2
    character(len=LENMEMPATH) :: pckMemPathsDummy
    real(DP) :: q
    !
    ! -- initialize counter
    idx = 0

    do i = 1, this%maxpackages
      ! -- Retrieve modelname1 and packagename1
      lloc = 1
      call urword(this%pckMemPaths(i), lloc, istart, istop, 1, ival, rval, -1, -1)
      pckMemPathsDummy = this%pckMemPaths(i)
      modelname1 = pckMemPathsDummy(istart:istop)
      call urword(this%pckMemPaths(i), lloc, istart, istop, 1, ival, rval, -1, -1)
      pckMemPathsDummy = this%pckMemPaths(i)
      packagename1 = pckMemPathsDummy(istart:istop)
      do j = 1, this%maxpackages
        ! -- Retrieve modelname2 and packagename2
        lloc = 1
        call urword(this%pckMemPaths(j), lloc, istart, istop, 1, ival, rval, &
                    -1, -1)
        pckMemPathsDummy = this%pckMemPaths(j)
        modelname2 = pckMemPathsDummy(istart:istop)
        call urword(this%pckMemPaths(j), lloc, istart, istop, 1, ival, rval, &
                    -1, -1)
        pckMemPathsDummy = this%pckMemPaths(j)
        packagename2 = pckMemPathsDummy(istart:istop)
        ipos = (i - 1) * this%maxpackages + j
        nitems = this%ientries(ipos)
        !
        ! -- nitems is the number of mover connections for this
        !    model-package / model-package combination.  Cycle if none.
        idx = idx + 1
        call this%budobj%budterm(idx)%reset(nitems)
        if (nitems < 1) cycle
        do n = 1, this%nmvr
          !
          ! -- pname1 is provider, pname2 is receiver
          !    flow is always negative because it is coming from provider
          if (this%pckMemPaths(i) == this%mvr(n)%mem_path_src) then
            if (this%pckMemPaths(j) == this%mvr(n)%mem_path_tgt) then
              !
              ! -- set q to qpactual
              q = -this%mvr(n)%qpactual
              !
              ! -- use mapped index (needed for lake to map outlet to lake number)
              n1 = this%mvr(n)%iRchNrSrcMapped
              !
              ! -- set receiver id to irch2
              n2 = this%mvr(n)%iRchNrTgt
              !
              ! -- check record into budget object
              call this%budobj%budterm(idx)%update_term(n1, n2, q)
            end if
          end if
        end do
      end do
    end do
    !
    ! --Terms are filled, now accumulate them for this time step
    call this%budobj%accumulate_terms()
  end subroutine fill_budobj

  !> @brief Set up output table
  !<
  subroutine mvr_setup_outputtab(this)
    ! -- dummy
    class(GwfMvrType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: ntabcol
    integer(I4B) :: ilen
    !
    ! -- allocate and initialize the output table
    if (this%iprflow /= 0) then
      !
      ! -- dimension table
      ntabcol = 7
      !
      ! -- initialize the output table object
      title = 'WATER MOVER PACKAGE ('//trim(this%packName)// &
              ') FLOW RATES'
      call table_cr(this%outputtab, this%packName, title)
      call this%outputtab%table_df(this%maxmvr, ntabcol, this%iout, &
                                   transient=.TRUE.)
      text = 'NUMBER'
      call this%outputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'PROVIDER LOCATION'
      ilen = LENMODELNAME + LENPACKAGENAME + 1
      call this%outputtab%initialize_column(text, ilen)
      text = 'PROVIDER ID'
      call this%outputtab%initialize_column(text, 10)
      text = 'AVAILABLE RATE'
      call this%outputtab%initialize_column(text, 10)
      text = 'PROVIDED RATE'
      call this%outputtab%initialize_column(text, 10)
      text = 'RECEIVER LOCATION'
      ilen = LENMODELNAME + LENPACKAGENAME + 1
      call this%outputtab%initialize_column(text, ilen)
      text = 'RECEIVER ID'
      call this%outputtab%initialize_column(text, 10)
      !
    end if
  end subroutine mvr_setup_outputtab

  !> @brief Set up output table
  !<
  subroutine mvr_print_outputtab(this)
    ! -- module
    use TdisModule, only: kstp, kper
    ! -- dummy
    class(GwfMvrType), intent(inout) :: this
    ! -- local
    character(len=LINELENGTH) :: title
    integer(I4B) :: i
    !
    ! -- set table kstp and kper
    call this%outputtab%set_kstpkper(kstp, kper)
    !
    ! -- Add terms and print the table
    title = 'WATER MOVER PACKAGE ('//trim(this%packName)// &
            ') FLOW RATES'
    call this%outputtab%set_title(title)
    call this%outputtab%set_maxbound(this%nmvr)
    do i = 1, this%nmvr
      call this%outputtab%add_term(i)
      call this%outputtab%add_term(this%mvr(i)%mem_path_src)
      call this%outputtab%add_term(this%mvr(i)%iRchNrSrc)
      call this%outputtab%add_term(this%mvr(i)%qavailable)
      call this%outputtab%add_term(this%mvr(i)%qpactual)
      call this%outputtab%add_term(this%mvr(i)%mem_path_tgt)
      call this%outputtab%add_term(this%mvr(i)%iRchNrTgt)
    end do
  end subroutine mvr_print_outputtab

  !> @brief Set mapped id
  !!
  !! For the budget output, we don't write outlet number,
  !! instead we write the lake number.  Normally the receiver
  !! number is the same as the feature number provided by the
  !! user.  For moving water from a lake, the user specifies the
  !! outlet number, not the lake number, in the mover package.
  !! The iRchNrSrcMapped variable contains the lake number, not
  !! the outlet number, and is written to the budget files.  For
  !! other packages, the iRchNrSrcMapped value is simply the well
  !! number, the stream reach, or the uzf cell number.
  !! This routine needs to be called each time a new set of movers
  !! is read.  It can't be called from within mvr_rp because the
  !! iprmap isn't updated by lake until lak_rp, which is called
  !! after mvr_rp.
  !<
  subroutine set_mapped_id(this)
    ! -- dummy
    class(GwfMvrType) :: this
    ! -- locals
    integer(I4B) :: i, mapped_id
    class(PackageMoverType), pointer :: pkg_mvr
    ! -- formats
    !
    ! -- set the feature maps
    allocate (pkg_mvr)
    do i = 1, this%nmvr
      call set_packagemover_pointer(pkg_mvr, this%mvr(i)%mem_path_src)
      mapped_id = pkg_mvr%iprmap(this%mvr(i)%iRchNrSrc)
      this%mvr(i)%iRchNrSrcMapped = mapped_id
    end do
    deallocate (pkg_mvr)
  end subroutine set_mapped_id

end module
