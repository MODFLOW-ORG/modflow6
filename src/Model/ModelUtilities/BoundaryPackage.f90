!> @brief This module contains the base boundary package
!!
!! This module contains the base model boundary package class that is
!! extended by all model boundary packages. The base model boundary
!! package extends the NumericalPackageType.
!<
module BndModule

  use KindModule, only: DP, LGP, I4B
  use ConstantsModule, only: LENAUXNAME, LENBOUNDNAME, LENFTYPE, &
                             DZERO, DONE, &
                             LENMODELNAME, LENPACKAGENAME, &
                             LENMEMPATH, MAXCHARLEN, LINELENGTH, &
                             DNODATA, LENLISTLABEL, LENPAKLOC, &
                             TABLEFT, TABCENTER
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, &
                       store_error_unit
  use NumericalPackageModule, only: NumericalPackageType
  use ObsModule, only: ObsType, obs_cr
  use TdisModule, only: delt, totimc
  use ObserveModule, only: ObserveType
  use InputOutputModule, only: GetUnit, openfile
  use TimeArraySeriesManagerModule, only: TimeArraySeriesManagerType
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType
  use ListModule, only: ListType
  use PackageMoverModule, only: PackageMoverType
  use BaseDisModule, only: DisBaseType
  use BlockParserModule, only: BlockParserType
  use TableModule, only: TableType, table_cr
  use CharacterStringModule, only: CharacterStringType
  use MatrixBaseModule

  implicit none

  private
  public :: BndType, AddBndToList, GetBndFromList
  public :: save_print_model_flows
  private :: CastAsBndClass

  !> @ brief BndType
  !!
  !!  Generic boundary package type.  This derived type can be overridden to
  !!  become concrete boundary package types.
  !<
  type, extends(NumericalPackageType) :: BndType
    ! -- characters
    character(len=LENLISTLABEL), pointer :: listlabel => null() !< title of table written for RP
    character(len=LENPACKAGENAME) :: text = '' !< text string for package flow term
    character(len=LENAUXNAME), dimension(:), pointer, &
      contiguous :: auxname => null() !< vector of auxname
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< copy of vector auxname that can be stored in memory manager
    character(len=LENBOUNDNAME), dimension(:), pointer, &
      contiguous :: boundname => null() !< vector of boundnames
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: boundname_cst => null() !< copy of vector boundname that can be stored in memory manager
    !
    ! -- scalars
    integer(I4B), pointer :: isadvpak => null() !< flag indicating package is advanced (1) or not (0)
    integer(I4B), pointer :: ibcnum => null() !< consecutive package number for this boundary condition
    integer(I4B), pointer :: maxbound => null() !< max number of boundaries
    integer(I4B), pointer :: nbound => null() !< number of boundaries for current stress period
    integer(I4B), pointer :: ncolbnd => null() !< number of columns of the bound array
    integer(I4B), pointer :: iscloc => null() !< bound column to scale with SFAC
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    integer(I4B), pointer :: inamedbound => null() !< flag to read boundnames
    integer(I4B), pointer :: iauxmultcol => null() !< column to use as multiplier for column iscloc
    integer(I4B), pointer :: npakeq => null() !< number of equations in this package (normally 0 unless package adds rows to matrix)
    integer(I4B), pointer :: ioffset => null() !< offset of this package in the model
    ! -- arrays
    integer(I4B), dimension(:), pointer, contiguous :: nodelist => null() !< vector of reduced node numbers
    integer(I4B), dimension(:), pointer, contiguous :: noupdateauxvar => null() !< override auxvars from being updated
    real(DP), dimension(:, :), pointer, contiguous :: bound => null() !< array of package specific boundary numbers
    real(DP), dimension(:), pointer, contiguous :: hcof => null() !< diagonal contribution
    real(DP), dimension(:), pointer, contiguous :: rhs => null() !< right-hand side contribution
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< auxiliary variable array
    real(DP), dimension(:), pointer, contiguous :: simvals => null() !< simulated values
    real(DP), dimension(:), pointer, contiguous :: simtomvr => null() !< simulated to mover values
    !
    ! -- water mover flag and object
    integer(I4B), pointer :: imover => null() !< flag indicating if the mover is active in the package
    type(PackageMoverType), pointer :: pakmvrobj => null() !< mover object for package
    !
    ! -- viscosity flag and safe-copy of conductance array
    integer(I4B), pointer :: ivsc => null() !< flag indicating if viscosity is active in the model
    real(DP), dimension(:), pointer, contiguous :: condinput => null() !< stores user-specified conductance values
    !
    ! -- timeseries
    type(TimeSeriesManagerType), pointer :: TsManager => null() !< time series manager
    type(TimeArraySeriesManagerType), pointer :: TasManager => null() !< time array series manager
    integer(I4B) :: indxconvertflux = 0 !< indxconvertflux is column of bound to multiply by area to convert flux to rate
    logical(LGP) :: AllowTimeArraySeries = .false.
    !
    ! -- pointers for observations
    integer(I4B), pointer :: inobspkg => null() !< unit number for obs package
    type(ObsType), pointer :: obs => null() !< observation package
    !
    ! -- pointers to model/solution variables
    integer(I4B), pointer :: neq !< number of equations for model
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< model ibound array
    real(DP), dimension(:), pointer, contiguous :: xnew => null() !< model dependent variable (head) for this time step
    real(DP), dimension(:), pointer, contiguous :: xold => null() !< model dependent variable for last time step
    real(DP), dimension(:), pointer, contiguous :: flowja => null() !< model intercell flows
    integer(I4B), dimension(:), pointer, contiguous :: icelltype => null() !< pointer to icelltype array in NPF
    character(len=LENMEMPATH) :: ictMemPath = '' !< memory path to the icelltype data (for GWF this is in NPF)
    !
    ! -- table objects
    type(TableType), pointer :: inputtab => null() !< input table object
    type(TableType), pointer :: outputtab => null() !< output table object for package flows writtent to the model listing file
    type(TableType), pointer :: errortab => null() !< package error table

  contains
    procedure, public :: bnd_df
    procedure :: bnd_ac
    procedure :: bnd_mc
    procedure :: bnd_ar
    procedure :: bnd_rp
    procedure :: bnd_ad
    procedure :: bnd_ck
    procedure :: bnd_reset
    procedure :: bnd_cf
    procedure :: bnd_fc
    procedure :: bnd_fn
    procedure :: bnd_nur
    procedure :: bnd_cc
    procedure :: bnd_cq
    procedure :: bnd_bd
    procedure :: bnd_ot_model_flows
    procedure :: bnd_ot_package_flows
    procedure :: bnd_ot_dv
    procedure :: bnd_ot_bdsummary
    procedure :: bnd_da

    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: pack_initialize
    procedure :: read_options => bnd_read_options
    procedure :: read_dimensions => bnd_read_dimensions
    procedure :: read_initial_attr => bnd_read_initial_attr
    procedure :: bnd_options
    procedure :: bnd_cq_simrate
    procedure :: bnd_cq_simtomvr
    procedure :: set_pointers
    procedure :: define_listlabel
    procedure :: copy_boundname
    procedure, private :: pak_setup_outputtab
    !
    ! -- procedures to support observations
    procedure, public :: bnd_obs_supported
    procedure, public :: bnd_df_obs
    procedure, public :: bnd_bd_obs
    procedure, public :: bnd_ot_obs
    procedure, public :: bnd_rp_obs
    !
    ! -- procedure to support time series
    procedure, public :: bnd_rp_ts
    !
    ! -- procedure to inform package that viscosity active
    procedure, public :: bnd_activate_viscosity
    !
    ! -- procedure to backup user-specified conductance
    procedure, private :: bnd_store_user_cond
    !
  end type BndType

contains

  !> @ brief Define boundary package options and dimensions
  !!
  !!  Define base boundary package options and dimensions for
  !!  a model boundary package.
  !<
  subroutine bnd_df(this, neq, dis)
    ! -- modules
    use TimeSeriesManagerModule, only: tsmanager_cr
    use TimeArraySeriesManagerModule, only: tasmanager_cr
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    integer(I4B), intent(inout) :: neq !< number of equations
    class(DisBaseType), pointer :: dis !< discretization object
    !
    ! -- set pointer to dis object for the model
    this%dis => dis
    !
    ! -- Create time series managers
    call tsmanager_cr(this%TsManager, this%iout)
    call tasmanager_cr(this%TasManager, dis, this%name_model, this%iout)
    !
    ! -- create obs package
    call obs_cr(this%obs, this%inobspkg)
    !
    ! -- Write information to model list file
    write (this%iout, 1) this%filtyp, trim(adjustl(this%text)), this%inunit
1   format(1X, /1X, a, ' -- ', a, ' PACKAGE, VERSION 8, 2/22/2014', &
           ' INPUT READ FROM UNIT ', I0)
    !
    ! -- Initialize block parser
    call this%parser%Initialize(this%inunit, this%iout)
    !
    ! -- set and read options
    call this%read_options()
    !
    ! -- Now that time series will have been read, need to call the df
    !    routine to define the manager
    call this%tsmanager%tsmanager_df()
    call this%tasmanager%tasmanager_df()
    !
    ! -- read the package dimensions block
    call this%read_dimensions()
    !
    ! -- update package moffset for packages that add rows
    if (this%npakeq > 0) then
      this%ioffset = neq - this%dis%nodes
    end if
    !
    ! -- update neq
    neq = neq + this%npakeq
    !
    ! -- Store information needed for observations
    if (this%bnd_obs_supported()) then
      call this%obs%obs_df(this%iout, this%packName, this%filtyp, this%dis)
      call this%bnd_df_obs()
    end if
  end subroutine bnd_df

  !> @ brief Add boundary package connection to matrix
  !!
  !!  Add boundary package connection to the matrix for packages that add
  !!  connections to the coefficient matrix. An example would be the GWF model
  !!  MAW package. Base implementation that must be extended.
  !<
  subroutine bnd_ac(this, moffset, sparse)
    ! -- modules
    use SparseModule, only: sparsematrix
    use SimModule, only: store_error
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    integer(I4B), intent(in) :: moffset !< solution matrix model offset
    type(sparsematrix), intent(inout) :: sparse !< sparse object
  end subroutine bnd_ac

  !> @ brief Map boundary package connection to matrix
  !!
  !!  Map boundary package connection to the matrix for packages that add
  !!  connections to the coefficient matrix. An example would be the GWF model
  !!  MAW package. Base implementation that must be extended.
  !<
  subroutine bnd_mc(this, moffset, matrix_sln)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    integer(I4B), intent(in) :: moffset !< solution matrix model offset
    class(MatrixBaseType), pointer :: matrix_sln !< global system matrix
  end subroutine bnd_mc

  !> @ brief Allocate and read method for boundary package
  !!
  !!  Generic method to allocate and read static data for model boundary
  !!  packages. A boundary package only needs to override this method if
  !!  input data varies from the standard boundary package.
  !<
  subroutine bnd_ar(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    !
    ! -- allocate and read observations
    call this%obs%obs_ar()
    !
    ! -- Allocate arrays in package superclass
    call this%allocate_arrays()
    !
    ! -- read optional initial package parameters
    call this%read_initial_attr()
    !
    ! -- setup pakmvrobj for standard stress packages
    if (this%imover == 1) then
      allocate (this%pakmvrobj)
      call this%pakmvrobj%ar(this%maxbound, 0, this%memoryPath)
    end if
  end subroutine bnd_ar

  !> @ brief Allocate and read method for package
  !!
  !!  Generic method to read and prepare period data for model boundary
  !!  packages. A boundary package only needs to override this method if
  !!  period data varies from the standard boundary package.
  !<
  subroutine bnd_rp(this)
    ! -- modules
    use TdisModule, only: kper, nper
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    integer(I4B) :: ierr
    integer(I4B) :: nlist
    logical(LGP) :: isfound
    character(len=LINELENGTH) :: line
    ! -- formats
    character(len=*), parameter :: fmtblkerr = &
      &"('Looking for BEGIN PERIOD iper.  Found ', a, ' instead.')"
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"
    character(len=*), parameter :: fmtnbd = &
      "(1X,/1X,'THE NUMBER OF ACTIVE ',A,'S (',I6, &
      &') IS GREATER THAN MAXIMUM(',I6,')')"
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
    ! -- read data if ionper == kper
    if (this%ionper == kper) then
      nlist = -1
      ! -- Remove all time-series and time-array-series links associated with
      !    this package.
      call this%TsManager%Reset(this%packName)
      call this%TasManager%Reset(this%packName)
      !
      ! -- Read data as a list
      call this%dis%read_list(this%parser%line_reader, &
                              this%parser%iuactive, this%iout, &
                              this%iprpak, nlist, this%inamedbound, &
                              this%iauxmultcol, this%nodelist, &
                              this%bound, this%auxvar, this%auxname, &
                              this%boundname, this%listlabel, &
                              this%packName, this%tsManager, this%iscloc)
      this%nbound = nlist
      !
      ! -- save user-specified conductance if vsc package is active
      if (this%ivsc == 1) then
        call this%bnd_store_user_cond(nlist, this%bound, this%condinput)
      end if
      !
      ! Define the tsLink%Text value(s) appropriately.
      ! E.g. for WEL package, entry 1, assign tsLink%Text = 'Q'
      ! For RIV package, entry 1 text = 'STAGE', entry 2 text = 'COND',
      !                  entry 3 text = 'RBOT'; etc.
      call this%bnd_rp_ts()
      !
      ! -- Terminate the block
      call this%parser%terminateblock()
      !
      ! -- Copy boundname into boundname_cst
      call this%copy_boundname()
      !
    else
      write (this%iout, fmtlsp) trim(this%filtyp)
    end if
  end subroutine bnd_rp

  !> @ brief Advance the boundary package
  !!
  !!  Advance data in the boundary package. The method sets advances
  !!  time series, time array series, and observation data. A boundary
  !!  package only needs to override this method if additional data
  !!  needs to be advanced.
  !<
  subroutine bnd_ad(this)
    ! -- dummy
    class(BndType) :: this !< BndType object
    ! -- local
    real(DP) :: begintime, endtime
    !
    ! -- Initialize time variables
    begintime = totimc
    endtime = begintime + delt
    !
    ! -- Advance the time series managers
    call this%TsManager%ad()
    call this%TasManager%ad()
    !
    ! -- For each observation, push simulated value and corresponding
    !    simulation time from "current" to "preceding" and reset
    !    "current" value.
    call this%obs%obs_ad()
  end subroutine bnd_ad

  !> @ brief Check boundary package period data
  !!
  !!  Check the boundary package period data. Base implementation that
  !!  must be extended by each model boundary package.
  !<
  subroutine bnd_ck(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    !
    ! -- check stress period data
    ! -- each package must override generic functionality
  end subroutine bnd_ck

  !> @ brief Reset bnd package before formulating
  !<
  subroutine bnd_reset(this)
    class(BndType) :: this !< BndType object

    if (this%imover == 1) then
      call this%pakmvrobj%reset()
    end if

  end subroutine bnd_reset

  !> @ brief Formulate the package hcof and rhs terms.
  !!
  !!  Formulate the hcof and rhs terms for the package that will be
  !!  added to the coefficient matrix and right-hand side vector.
  !!  Base implementation that must be extended by each model
  !!  boundary package.
  !<
  subroutine bnd_cf(this)
    ! -- modules
    class(BndType) :: this !< BndType object
    !
    ! -- bnd has no cf routine
  end subroutine bnd_cf

  !> @ brief Copy hcof and rhs terms into solution.
  !!
  !!  Add the hcof and rhs terms for the boundary package to the
  !!  coefficient matrix and right-hand side vector. A boundary
  !!  package only needs to override this method if it is different for
  !!  a specific boundary package.
  !<
  subroutine bnd_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(BndType) :: this !< BndType object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n
    integer(I4B) :: ipos
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
    end do
  end subroutine bnd_fc

  !> @ brief Add Newton-Raphson terms for package into solution.
  !!
  !!  Calculate and add the Newton-Raphson terms for the boundary package
  !!  to the coefficient matrix and right-hand side vector. A boundary
  !!  package only needs to override this method if a specific boundary
  !!  package needs to add Newton-Raphson terms.
  !<
  subroutine bnd_fn(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(BndType) :: this !< BndType object
    real(DP), dimension(:), intent(inout) :: rhs !< right-hand side vector for model
    integer(I4B), dimension(:), intent(in) :: ia !< solution CRS row pointers
    integer(I4B), dimension(:), intent(in) :: idxglo !< mapping vector for model (local) to solution (global)
    class(MatrixBaseType), pointer :: matrix_sln !< solution coefficient matrix
    !
    ! -- No addition terms for Newton-Raphson with constant conductance
    !    boundary conditions
  end subroutine bnd_fn

  !> @ brief Apply Newton-Raphson under-relaxation for package.
  !!
  !!  Apply Newton-Raphson under-relaxation for a boundary package. A boundary
  !!  package only needs to override this method if a specific boundary
  !!  package needs to apply Newton-Raphson under-relaxation. An example is
  !!  the MAW package which adds rows to the system of equations and may need
  !!  to have the dependent-variable constrained by the bottom of the model.
  !<
  subroutine bnd_nur(this, neqpak, x, xtemp, dx, inewtonur, dxmax, locmax)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    integer(I4B), intent(in) :: neqpak !< number of equations in the package
    real(DP), dimension(neqpak), intent(inout) :: x !< dependent variable
    real(DP), dimension(neqpak), intent(in) :: xtemp !< previous dependent variable
    real(DP), dimension(neqpak), intent(inout) :: dx !< change in dependent variable
    integer(I4B), intent(inout) :: inewtonur !< flag indicating if newton-raphson under-relaxation should be applied
    real(DP), intent(inout) :: dxmax !< maximum change in the dependent variable
    integer(I4B), intent(inout) :: locmax !< location of the maximum change in the dependent variable
    ! -- local
    !
    ! -- Newton-Raphson under-relaxation
  end subroutine bnd_nur

  !> @ brief Convergence check for package.
  !!
  !!  Perform additional convergence checks on the flow between the package
  !!  and the model it is attached to. This additional convergence check is
  !!  applied to packages that solve their own continuity equation as
  !!  part of the formulate step at the beginning of a Picard iteration.
  !!  A boundary package only needs to override this method if a specific boundary
  !!  package solves its own continuity equation. Example packages that implement
  !!  this additional convergence check is the CSUB, SFR, LAK, and UZF packages.
  !<
  subroutine bnd_cc(this, innertot, kiter, iend, icnvgmod, cpak, ipak, dpak)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    integer(I4B), intent(in) :: innertot !< total number of inner iterations
    integer(I4B), intent(in) :: kiter !< Picard iteration number
    integer(I4B), intent(in) :: iend !< flag indicating if this is the last Picard iteration
    integer(I4B), intent(in) :: icnvgmod !< flag inficating if the model has met specific convergence criteria
    character(len=LENPAKLOC), intent(inout) :: cpak !< string for user node
    integer(I4B), intent(inout) :: ipak !< location of the maximum dependent variable change
    real(DP), intent(inout) :: dpak !< maximum dependent variable change
    !
    ! -- No addition convergence check for boundary conditions
  end subroutine bnd_cc

  !> @ brief Calculate advanced package flows.
  !!
  !!  Calculate the flow between connected advanced package control volumes.
  !!  Only advanced boundary packages need to override this method.
  !<
  subroutine bnd_cq(this, x, flowja, iadv)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    real(DP), dimension(:), intent(in) :: x !< current dependent-variable value
    real(DP), dimension(:), contiguous, intent(inout) :: flowja !< flow between two connected control volumes
    integer(I4B), optional, intent(in) :: iadv !< flag that indicates if this is an advance package
    ! -- local
    integer(I4B) :: imover
    !
    ! -- check for iadv optional variable to indicate this is an advanced
    !    package and that mover calculations should not be done here
    if (present(iadv)) then
      if (iadv == 1) then
        imover = 0
      else
        imover = 1
      end if
    else
      imover = this%imover
    end if
    !
    ! -- Calculate package flows.  In the first call, simval is calculated
    !    from hcof, rhs, and head.  The second call may reduce the value in
    !    simval by what is sent to the mover.  The mover rate is stored in
    !    simtomvr.  imover is set to zero here for advanced packages, which
    !    handle and store mover quantities separately.
    call this%bnd_cq_simrate(x, flowja, imover)
    if (imover == 1) then
      call this%bnd_cq_simtomvr(flowja)
    end if
  end subroutine bnd_cq

  !> @ brief Calculate simrate.
  !!
  !!  Calculate the flow between package and the model (for example, GHB and
  !!  groundwater cell) and store in the simvals variable. This method only
  !!  needs to be overridden if a different calculation needs to be made.
  !<
  subroutine bnd_cq_simrate(this, hnew, flowja, imover)
    ! -- dummy
    class(BndType) :: this !< BndType object
    real(DP), dimension(:), intent(in) :: hnew !< current dependent-variable value
    real(DP), dimension(:), intent(inout) :: flowja !< flow between package and model
    integer(I4B), intent(in) :: imover !< flag indicating if the mover package is active
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: idiag
    real(DP) :: rrate
    !
    ! -- If no boundaries, skip flow calculations.
    if (this%nbound > 0) then
      !
      ! -- Loop through each boundary calculating flow.
      do i = 1, this%nbound
        node = this%nodelist(i)
        !
        ! -- If cell is no-flow or constant-head, then ignore it.
        rrate = DZERO
        if (node > 0) then
          idiag = this%dis%con%ia(node)
          if (this%ibound(node) > 0) then
            !
            ! -- Calculate the flow rate into the cell.
            rrate = this%hcof(i) * hnew(node) - this%rhs(i)
          end if
          flowja(idiag) = flowja(idiag) + rrate
        end if
        !
        ! -- Save simulated value to simvals array.
        this%simvals(i) = rrate
        !
      end do
    end if
  end subroutine bnd_cq_simrate

  !> @ brief Calculate flow to the mover.
  !!
  !!  Calculate the flow between package and the model that is sent to the
  !!  mover package and store in the simtomvr variable. This method only
  !!  needs to be overridden if a different calculation needs to be made.
  !<
  subroutine bnd_cq_simtomvr(this, flowja)
    ! -- dummy
    class(BndType) :: this !< BndType object
    real(DP), dimension(:), intent(inout) :: flowja !< flow between package and model
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: node
    real(DP) :: q
    real(DP) :: fact
    real(DP) :: rrate
    !
    ! -- If no boundaries, skip flow calculations.
    if (this%nbound > 0) then
      !
      ! -- Loop through each boundary calculating flow.
      do i = 1, this%nbound
        node = this%nodelist(i)
        !
        ! -- If cell is no-flow or constant-head, then ignore it.
        rrate = DZERO
        if (node > 0) then
          if (this%ibound(node) > 0) then
            !
            ! -- Calculate the flow rate into the cell.
            q = this%simvals(i)

            if (q < DZERO) then
              rrate = this%pakmvrobj%get_qtomvr(i)
              !
              ! -- Evaluate if qtomvr exceeds the calculated rrate.
              !    When fact is greater than 1, qtomvr is numerically
              !    larger than rrate (which should never happen) and
              !    represents a water budget error. When this happens,
              !    rrate is set to 0. so that the water budget error is
              !    correctly accounted for in the listing water budget.
              fact = -rrate / q
              if (fact > DONE) then
                ! -- all flow goes to mover
                q = DZERO
              else
                ! -- magnitude of rrate (which is negative) is reduced by
                !    qtomvr (which is positive)
                q = q + rrate
              end if
              this%simvals(i) = q

              if (rrate > DZERO) then
                rrate = -rrate
              end if
            end if
          end if
        end if
        !
        ! -- Save simulated value to simtomvr array.
        this%simtomvr(i) = rrate
        !
      end do
    end if
  end subroutine bnd_cq_simtomvr

  !> @ brief Add package flows to model budget.
  !!
  !!  Add the flow between package and the model (ratin and ratout) to the
  !!  model budget. This method only needs to be overridden if a different
  !!  calculation needs to be made.
  !<
  subroutine bnd_bd(this, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(BndType) :: this !< BndType object
    type(BudgetType), intent(inout) :: model_budget !< model budget object
    ! -- local
    character(len=LENPACKAGENAME) :: text
    real(DP) :: ratin
    real(DP) :: ratout
    integer(I4B) :: isuppress_output
    !
    ! -- initialize local variables
    isuppress_output = 0
    !
    ! -- call accumulator and add to the model budget
    call rate_accumulator(this%simvals(1:this%nbound), ratin, ratout)
    call model_budget%addentry(ratin, ratout, delt, this%text, &
                               isuppress_output, this%packName)
    if (this%imover == 1 .and. this%isadvpak == 0) then
      text = trim(adjustl(this%text))//'-TO-MVR'
      text = adjustr(text)
      call rate_accumulator(this%simtomvr(1:this%nbound), ratin, ratout)
      call model_budget%addentry(ratin, ratout, delt, text, &
                                 isuppress_output, this%packName)
    end if
  end subroutine bnd_bd

  !> @ brief Output advanced package flow terms.
  !!
  !!  Output advanced boundary package flow terms. This method only needs to
  !!  be overridden for advanced packages that save flow terms than contribute
  !!  to the continuity equation for each control volume.
  !<
  subroutine bnd_ot_package_flows(this, icbcfl, ibudfl)
    ! -- dummy
    class(BndType) :: this !< BndType object
    integer(I4B), intent(in) :: icbcfl !< flag and unit number for cell-by-cell output
    integer(I4B), intent(in) :: ibudfl !< flag indication if cell-by-cell data should be saved
    !
    ! -- override for advanced packages
  end subroutine bnd_ot_package_flows

  !> @ brief Output advanced package dependent-variable terms.
  !!
  !!  Output advanced boundary package dependent-variable terms. This method only needs
  !!  to be overridden for advanced packages that save dependent variable terms
  !!  for each control volume.
  !<
  subroutine bnd_ot_dv(this, idvsave, idvprint)
    ! -- dummy
    class(BndType) :: this !< BndType object
    integer(I4B), intent(in) :: idvsave !< flag and unit number for dependent-variable output
    integer(I4B), intent(in) :: idvprint !< flag indicating if dependent-variable should be written to the model listing file
    !
    ! -- override for advanced packages
  end subroutine bnd_ot_dv

  !> @ brief Output advanced package budget summary.
  !!
  !!  Output advanced boundary package budget summary. This method only needs
  !!  to be overridden for advanced packages that save budget summaries
  !!  to the model listing file.
  !<
  subroutine bnd_ot_bdsummary(this, kstp, kper, iout, ibudfl)
    ! -- dummy
    class(BndType) :: this !< BndType object
    integer(I4B), intent(in) :: kstp !< time step number
    integer(I4B), intent(in) :: kper !< period number
    integer(I4B), intent(in) :: iout !< flag and unit number for the model listing file
    integer(I4B), intent(in) :: ibudfl !< flag indicating budget should be written
    !
    ! -- override for advanced packages
  end subroutine bnd_ot_bdsummary

  !> @ brief Output package flow terms.
  !!
  !!  Output flow terms between the boundary package and model to a binary file and/or
  !!  print flows to the model listing file. This method should not need to
  !!  be overridden.
  !<
  subroutine bnd_ot_model_flows(this, icbcfl, ibudfl, icbcun, imap)
    ! -- dummy
    class(BndType) :: this !< BndType object
    integer(I4B), intent(in) :: icbcfl !< flag for cell-by-cell output
    integer(I4B), intent(in) :: ibudfl !< flag indication if cell-by-cell data should be saved
    integer(I4B), intent(in) :: icbcun !< unit number for cell-by-cell output
    integer(I4B), dimension(:), optional, intent(in) :: imap !< mapping vector that converts the 1 to nbound values to lake number, maw number, etc.
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LENPACKAGENAME) :: text
    integer(I4B) :: imover
    !
    ! -- Call generic subroutine to save and print simvals and simtomvr
    title = trim(adjustl(this%text))//' PACKAGE ('//trim(this%packName)// &
            ') FLOW RATES'
    if (present(imap)) then
      call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                  this%outputtab, this%nbound, this%nodelist, &
                                  this%simvals, this%ibound, title, this%text, &
                                  this%ipakcb, this%dis, this%naux, &
                                  this%name_model, this%name_model, &
                                  this%name_model, this%packName, &
                                  this%auxname, this%auxvar, this%iout, &
                                  this%inamedbound, this%boundname, imap)
    else
      call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                  this%outputtab, this%nbound, this%nodelist, &
                                  this%simvals, this%ibound, title, this%text, &
                                  this%ipakcb, this%dis, this%naux, &
                                  this%name_model, this%name_model, &
                                  this%name_model, this%packName, &
                                  this%auxname, this%auxvar, this%iout, &
                                  this%inamedbound, this%boundname)
    end if
    !
    ! -- Set mover flag, and shut off if this is an advanced package.  Advanced
    !    packages must handle mover flows differently by including them in
    !    their balance equations.  These simtomvr flows are the general
    !    flow to mover terms calculated by bnd_cq_simtomvr()
    imover = this%imover
    if (this%isadvpak /= 0) imover = 0
    if (imover == 1) then
      text = trim(adjustl(this%text))//'-TO-MVR'
      text = adjustr(text)
      title = trim(adjustl(this%text))//' PACKAGE ('// &
              trim(this%packName)//') FLOW RATES TO-MVR'
      call save_print_model_flows(icbcfl, ibudfl, icbcun, this%iprflow, &
                                  this%outputtab, this%nbound, this%nodelist, &
                                  this%simtomvr, this%ibound, title, text, &
                                  this%ipakcb, this%dis, this%naux, &
                                  this%name_model, this%name_model, &
                                  this%name_model, this%packName, &
                                  this%auxname, this%auxvar, this%iout, &
                                  this%inamedbound, this%boundname)
    end if
  end subroutine bnd_ot_model_flows

  !> @ brief Deallocate package memory
  !!
  !!  Deallocate base boundary package scalars and arrays. This method
  !!  only needs to be overridden if additional variables are defined
  !!  for a specific package.
  !<
  subroutine bnd_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(BndType) :: this !< BndType object
    !
    ! -- deallocate arrays
    call mem_deallocate(this%nodelist, 'NODELIST', this%memoryPath)
    call mem_deallocate(this%noupdateauxvar, 'NOUPDATEAUXVAR', this%memoryPath)
    call mem_deallocate(this%bound, 'BOUND', this%memoryPath)
    call mem_deallocate(this%condinput, 'CONDINPUT', this%memoryPath)
    call mem_deallocate(this%hcof, 'HCOF', this%memoryPath)
    call mem_deallocate(this%rhs, 'RHS', this%memoryPath)
    call mem_deallocate(this%simvals, 'SIMVALS', this%memoryPath)
    call mem_deallocate(this%simtomvr, 'SIMTOMVR', this%memoryPath)
    call mem_deallocate(this%auxvar, 'AUXVAR', this%memoryPath)
    call mem_deallocate(this%boundname, 'BOUNDNAME', this%memoryPath)
    call mem_deallocate(this%boundname_cst, 'BOUNDNAME_CST', this%memoryPath)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxname_cst, 'AUXNAME_CST', this%memoryPath)
    nullify (this%icelltype)
    !
    ! -- pakmvrobj
    if (this%imover /= 0) then
      call this%pakmvrobj%da()
      deallocate (this%pakmvrobj)
      nullify (this%pakmvrobj)
    end if
    !
    ! -- input table object
    if (associated(this%inputtab)) then
      call this%inputtab%table_da()
      deallocate (this%inputtab)
      nullify (this%inputtab)
    end if
    !
    ! -- output table object
    if (associated(this%outputtab)) then
      call this%outputtab%table_da()
      deallocate (this%outputtab)
      nullify (this%outputtab)
    end if
    !
    ! -- error table object
    if (associated(this%errortab)) then
      call this%errortab%table_da()
      deallocate (this%errortab)
      nullify (this%errortab)
    end if
    !
    ! -- deallocate character variables
    call mem_deallocate(this%listlabel, 'LISTLABEL', this%memoryPath)
    !
    ! -- Deallocate scalars
    call mem_deallocate(this%isadvpak)
    call mem_deallocate(this%ibcnum)
    call mem_deallocate(this%maxbound)
    call mem_deallocate(this%nbound)
    call mem_deallocate(this%ncolbnd)
    call mem_deallocate(this%iscloc)
    call mem_deallocate(this%naux)
    call mem_deallocate(this%inamedbound)
    call mem_deallocate(this%iauxmultcol)
    call mem_deallocate(this%inobspkg)
    call mem_deallocate(this%imover)
    call mem_deallocate(this%npakeq)
    call mem_deallocate(this%ioffset)
    call mem_deallocate(this%ivsc)
    !
    ! -- deallocate methods on objects
    call this%obs%obs_da()
    call this%TsManager%da()
    call this%TasManager%da()
    !
    ! -- deallocate objects
    deallocate (this%obs)
    deallocate (this%TsManager)
    deallocate (this%TasManager)
    nullify (this%TsManager)
    nullify (this%TasManager)
    !
    ! -- Deallocate parent object
    call this%NumericalPackageType%da()
  end subroutine bnd_da

  !> @ brief Allocate package scalars
  !!
  !!  Allocate and initialize base boundary package scalars. This method
  !!  only needs to be overridden if additional scalars are defined
  !!  for a specific package.
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(BndType) :: this !< BndType object
    ! -- local
    integer(I4B), pointer :: imodelnewton => null()
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- allocate character variables
    call mem_allocate(this%listlabel, LENLISTLABEL, 'LISTLABEL', &
                      this%memoryPath)
    !
    ! -- allocate integer variables
    call mem_allocate(this%isadvpak, 'ISADVPAK', this%memoryPath)
    call mem_allocate(this%ibcnum, 'IBCNUM', this%memoryPath)
    call mem_allocate(this%maxbound, 'MAXBOUND', this%memoryPath)
    call mem_allocate(this%nbound, 'NBOUND', this%memoryPath)
    call mem_allocate(this%ncolbnd, 'NCOLBND', this%memoryPath)
    call mem_allocate(this%iscloc, 'ISCLOC', this%memoryPath)
    call mem_allocate(this%naux, 'NAUX', this%memoryPath)
    call mem_allocate(this%inamedbound, 'INAMEDBOUND', this%memoryPath)
    call mem_allocate(this%iauxmultcol, 'IAUXMULTCOL', this%memoryPath)
    call mem_allocate(this%inobspkg, 'INOBSPKG', this%memoryPath)
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%imover, 'IMOVER', this%memoryPath)
    !
    ! -- allocate flag for determining if vsc active
    call mem_allocate(this%ivsc, 'IVSC', this%memoryPath)
    !
    ! -- allocate scalars for packages that add rows to the matrix (e.g. MAW)
    call mem_allocate(this%npakeq, 'NPAKEQ', this%memoryPath)
    call mem_allocate(this%ioffset, 'IOFFSET', this%memoryPath)
    !
    ! -- allocate TS objects
    allocate (this%TsManager)
    allocate (this%TasManager)
    !
    ! -- allocate text strings
    call mem_allocate(this%auxname, LENAUXNAME, 0, 'AUXNAME', this%memoryPath)
    call mem_allocate(this%auxname_cst, LENAUXNAME, 0, 'AUXNAME_CST', &
                      this%memoryPath)
    !
    ! -- Initialize variables
    this%isadvpak = 0
    this%ibcnum = 0
    this%maxbound = 0
    this%nbound = 0
    this%ncolbnd = 0
    this%iscloc = 0
    this%naux = 0
    this%inamedbound = 0
    this%iauxmultcol = 0
    this%inobspkg = 0
    this%imover = 0
    this%npakeq = 0
    this%ioffset = 0
    this%ivsc = 0
    !
    ! -- Set pointer to model inewton variable
    call mem_setptr(imodelnewton, 'INEWTON', create_mem_path(this%name_model))
    this%inewton = imodelnewton
    imodelnewton => null()
  end subroutine allocate_scalars

  !> @ brief Allocate package arrays
  !!
  !!  Allocate and initialize base boundary package arrays. This method
  !!  only needs to be overridden if additional arrays are defined
  !!  for a specific package.
  !<
  subroutine allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(BndType) :: this !< BndType object
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist !< package nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar !< package aux variable array
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    !
    ! -- Point nodelist if it is passed in, otherwise allocate
    if (present(nodelist)) then
      this%nodelist => nodelist
    else
      call mem_allocate(this%nodelist, this%maxbound, 'NODELIST', &
                        this%memoryPath)
      do j = 1, this%maxbound
        this%nodelist(j) = 0
      end do
    end if
    !
    ! -- noupdateauxvar (allows an external caller to stop auxvars from being
    !    recalculated
    call mem_allocate(this%noupdateauxvar, this%naux, 'NOUPDATEAUXVAR', &
                      this%memoryPath)
    this%noupdateauxvar(:) = 0
    !
    ! -- Allocate the bound array
    call mem_allocate(this%bound, this%ncolbnd, this%maxbound, 'BOUND', &
                      this%memoryPath)
    !
    !-- Allocate array for storing user-specified conductances
    !   Will be reallocated to size maxbound if vsc active
    call mem_allocate(this%condinput, 0, 'CONDINPUT', this%memoryPath)
    !
    ! -- Allocate hcof and rhs
    call mem_allocate(this%hcof, this%maxbound, 'HCOF', this%memoryPath)
    call mem_allocate(this%rhs, this%maxbound, 'RHS', this%memoryPath)
    !
    ! -- Allocate the simvals array
    call mem_allocate(this%simvals, this%maxbound, 'SIMVALS', this%memoryPath)
    if (this%imover == 1) then
      call mem_allocate(this%simtomvr, this%maxbound, 'SIMTOMVR', &
                        this%memoryPath)
      do i = 1, this%maxbound
        this%simtomvr(i) = DZERO
      end do
    else
      call mem_allocate(this%simtomvr, 0, 'SIMTOMVR', this%memoryPath)
    end if
    !
    ! -- Point or allocate auxvar
    if (present(auxvar)) then
      this%auxvar => auxvar
    else
      call mem_allocate(this%auxvar, this%naux, this%maxbound, 'AUXVAR', &
                        this%memoryPath)
      do i = 1, this%maxbound
        do j = 1, this%naux
          this%auxvar(j, i) = DZERO
        end do
      end do
    end if
    !
    ! -- Allocate boundname
    if (this%inamedbound /= 0) then
      call mem_allocate(this%boundname, LENBOUNDNAME, this%maxbound, &
                        'BOUNDNAME', this%memoryPath)
      call mem_allocate(this%boundname_cst, LENBOUNDNAME, this%maxbound, &
                        'BOUNDNAME_CST', this%memoryPath)
    else
      call mem_allocate(this%boundname, LENBOUNDNAME, 0, &
                        'BOUNDNAME', this%memoryPath)
      call mem_allocate(this%boundname_cst, LENBOUNDNAME, 0, &
                        'BOUNDNAME_CST', this%memoryPath)
    end if
    !
    ! -- Set pointer to ICELLTYPE. For GWF boundary packages,
    !    this%ictMemPath will be 'NPF'.  If boundary packages do not set
    !    this%ictMemPath, then icelltype will remain as null()
    if (this%ictMemPath /= '') then
      call mem_setptr(this%icelltype, 'ICELLTYPE', this%ictMemPath)
    end if
    !
    ! -- Initialize values
    do j = 1, this%maxbound
      do i = 1, this%ncolbnd
        this%bound(i, j) = DZERO
      end do
    end do
    do i = 1, this%maxbound
      this%hcof(i) = DZERO
      this%rhs(i) = DZERO
    end do
    !
    ! -- setup the output table
    call this%pak_setup_outputtab()
  end subroutine allocate_arrays

  !> @ brief Allocate and initialize select package members
  !!
  !!  Allocate and initialize select base boundary package members.
  !!  This method needs to be overridden by a package if it is
  !!  needed for a specific package.
  !<
  subroutine pack_initialize(this)
    ! -- dummy
    class(BndType) :: this !< BndType object
  end subroutine pack_initialize

  !> @ brief Set pointers to model variables
  !!
  !!  Set pointers to model variables so that a package has access to these
  !!  variables. This base method should not need to be overridden.
  !<
  subroutine set_pointers(this, neq, ibound, xnew, xold, flowja)
    ! -- dummy
    class(BndType) :: this !< BndType object
    integer(I4B), pointer :: neq !< number of equations in the model
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< model idomain
    real(DP), dimension(:), pointer, contiguous :: xnew !< current dependent variable
    real(DP), dimension(:), pointer, contiguous :: xold !< previous dependent variable
    real(DP), dimension(:), pointer, contiguous :: flowja !< connection flow terms
    !
    ! -- Set the pointers
    this%neq => neq
    this%ibound => ibound
    this%xnew => xnew
    this%xold => xold
    this%flowja => flowja
  end subroutine set_pointers

  !> @ brief Read additional options for package
  !!
  !!  Read base options for boundary packages.
  !<
  subroutine bnd_read_options(this)
    ! -- modules
    use InputOutputModule, only: urdaux
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    character(len=:), allocatable :: line
    character(len=LINELENGTH) :: fname
    character(len=LINELENGTH) :: keyword
    character(len=LENAUXNAME) :: sfacauxname
    character(len=LENAUXNAME), dimension(:), allocatable :: caux
    integer(I4B) :: lloc
    integer(I4B) :: istart
    integer(I4B) :: istop
    integer(I4B) :: n
    integer(I4B) :: ierr
    integer(I4B) :: inobs
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    logical(LGP) :: foundchildclassoption
    ! -- format
    character(len=*), parameter :: fmtflow = &
      &"(4x, 'FLOWS WILL BE SAVED TO FILE: ', a, /4x, 'OPENED ON UNIT: ', I7)"
    character(len=*), parameter :: fmtflow2 = &
      &"(4x, 'FLOWS WILL BE SAVED TO BUDGET FILE SPECIFIED IN OUTPUT CONTROL')"
    character(len=*), parameter :: fmttas = &
      &"(4x, 'TIME-ARRAY SERIES DATA WILL BE READ FROM FILE: ', a)"
    character(len=*), parameter :: fmtts = &
      &"(4x, 'TIME-SERIES DATA WILL BE READ FROM FILE: ', a)"
    character(len=*), parameter :: fmtnme = &
      &"(a, i0, a)"
    !
    ! -- set default options
    !
    ! -- get options block
    call this%parser%GetBlock('OPTIONS', isfound, ierr, &
                              supportOpenClose=.true., blockRequired=.false.)
    !
    ! -- parse options block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
        //' OPTIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) then
          exit
        end if
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('AUX', 'AUXILIARY')
          call this%parser%GetRemainingLine(line)
          lloc = 1
          call urdaux(this%naux, this%parser%iuactive, this%iout, lloc, &
                      istart, istop, caux, line, this%text)
          call mem_reallocate(this%auxname, LENAUXNAME, this%naux, &
                              'AUXNAME', this%memoryPath)
          call mem_reallocate(this%auxname_cst, LENAUXNAME, this%naux, &
                              'AUXNAME_CST', this%memoryPath)
          do n = 1, this%naux
            this%auxname(n) = caux(n)
            this%auxname_cst(n) = caux(n)
          end do
          deallocate (caux)
        case ('SAVE_FLOWS')
          this%ipakcb = -1
          write (this%iout, fmtflow2)
        case ('PRINT_INPUT')
          this%iprpak = 1
          write (this%iout, '(4x,a)') &
            'LISTS OF '//trim(adjustl(this%text))//' CELLS WILL BE PRINTED.'
        case ('PRINT_FLOWS')
          this%iprflow = 1
          write (this%iout, '(4x,a)') trim(adjustl(this%text))// &
            ' FLOWS WILL BE PRINTED TO LISTING FILE.'
        case ('BOUNDNAMES')
          this%inamedbound = 1
          write (this%iout, '(4x,a)') trim(adjustl(this%text))// &
            ' BOUNDARIES HAVE NAMES IN LAST COLUMN.'
        case ('TS6')
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'TS6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
          end if
          call this%parser%GetString(fname)
          write (this%iout, fmtts) trim(fname)
          call this%TsManager%add_tsfile(fname, this%inunit)
        case ('TAS6')
          if (this%AllowTimeArraySeries) then
            if (.not. this%dis%supports_layers()) then
              errmsg = 'TAS6 FILE cannot be used '// &
                       'with selected discretization type.'
              call store_error(errmsg)
            end if
          else
            errmsg = 'The '//trim(this%filtyp)// &
                     ' package does not support TIMEARRAYSERIESFILE'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'TAS6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          call this%parser%GetString(fname)
          write (this%iout, fmttas) trim(fname)
          call this%TasManager%add_tasfile(fname)
        case ('AUXMULTNAME')
          call this%parser%GetStringCaps(sfacauxname)
          this%iauxmultcol = -1
          write (this%iout, '(4x,a,a)') &
            'AUXILIARY MULTIPLIER NAME: ', sfacauxname
        case ('OBS6')
          call this%parser%GetStringCaps(keyword)
          if (trim(adjustl(keyword)) /= 'FILEIN') then
            errmsg = 'OBS6 keyword must be followed by "FILEIN" '// &
                     'then by filename.'
            call store_error(errmsg)
            call this%parser%StoreErrorUnit()
          end if
          if (this%obs%active) then
            errmsg = 'Multiple OBS6 keywords detected in OPTIONS block. '// &
                     'Only one OBS6 entry allowed for a package.'
            call store_error(errmsg)
          end if
          this%obs%active = .true.
          call this%parser%GetString(this%obs%inputFilename)
          inobs = GetUnit()
          call openfile(inobs, this%iout, this%obs%inputFilename, 'OBS')
          this%obs%inUnitObs = inobs
          !
          ! -- right now these are options that are only available in the
          !    development version and are not included in the documentation.
          !    These options are only available when IDEVELOPMODE in
          !    constants module is set to 1
        case ('DEV_NO_NEWTON')
          call this%parser%DevOpt()
          this%inewton = 0
          write (this%iout, '(4x,a)') &
            'NEWTON-RAPHSON method disabled for unconfined cells'
        case default
          !
          ! -- Check for child class options
          call this%bnd_options(keyword, foundchildclassoption)
          !
          ! -- No child class options found, so print error message
          if (.not. foundchildclassoption) then
            write (errmsg, '(a,3(1x,a))') &
              'UNKNOWN', trim(adjustl(this%text)), 'OPTION:', trim(keyword)
            call store_error(errmsg)
          end if
        end select
      end do
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' OPTIONS'
    else
      write (this%iout, '(1x,a)') 'NO '//trim(adjustl(this%text))// &
        ' OPTION BLOCK DETECTED.'
    end if
    !
    ! -- AUXMULTNAME was specified, so find column of auxvar that will be multiplier
    if (this%iauxmultcol < 0) then
      !
      ! -- Error if no aux variable specified
      if (this%naux == 0) then
        write (errmsg, '(a,2(1x,a))') &
          'AUXMULTNAME was specified as', trim(adjustl(sfacauxname)), &
          'but no AUX variables specified.'
        call store_error(errmsg)
      end if
      !
      ! -- Assign mult column
      this%iauxmultcol = 0
      do n = 1, this%naux
        if (sfacauxname == this%auxname(n)) then
          this%iauxmultcol = n
          exit
        end if
      end do
      !
      ! -- Error if aux variable cannot be found
      if (this%iauxmultcol == 0) then
        write (errmsg, '(a,2(1x,a))') &
          'AUXMULTNAME was specified as', trim(adjustl(sfacauxname)), &
          'but no AUX variable found with this name.'
        call store_error(errmsg)
      end if
    end if
    !
    ! -- terminate if errors were detected
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end subroutine bnd_read_options

  !> @ brief Read dimensions for package
  !!
  !!  Read base dimensions for boundary packages. This method should not
  !!  need to be overridden unless more than MAXBOUND is specified in the
  !!  DIMENSIONS block.
  !<
  subroutine bnd_read_dimensions(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    character(len=LINELENGTH) :: keyword
    logical(LGP) :: isfound
    logical(LGP) :: endOfBlock
    integer(I4B) :: ierr
    !
    ! -- get dimensions block
    call this%parser%GetBlock('DIMENSIONS', isfound, ierr, &
                              supportOpenClose=.true.)
    !
    ! -- parse dimensions block if detected
    if (isfound) then
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text))// &
        ' DIMENSIONS'
      do
        call this%parser%GetNextLine(endOfBlock)
        if (endOfBlock) exit
        call this%parser%GetStringCaps(keyword)
        select case (keyword)
        case ('MAXBOUND')
          this%maxbound = this%parser%GetInteger()
          write (this%iout, '(4x,a,i7)') 'MAXBOUND = ', this%maxbound
        case default
          write (errmsg, '(a,3(1x,a))') &
            'Unknown', trim(this%text), 'dimension:', trim(keyword)
          call store_error(errmsg)
        end select
      end do
      !
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
    else
      call store_error('Required DIMENSIONS block not found.')
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- verify dimensions were set
    if (this%maxbound <= 0) then
      write (errmsg, '(a)') 'MAXBOUND must be an integer greater than zero.'
      call store_error(errmsg)
    end if
    !
    ! -- terminate if there are errors
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
  end subroutine bnd_read_dimensions

  !> @ brief Store user-specified conductances when vsc is active
  !!
  !!  VSC will update boundary package conductance values.  Because
  !!  viscosity can change every stress period, but user-specified
  !!  conductances may not, the base user-input should be stored in
  !!  backup array so that viscosity-updated conductances may be
  !!  recalculated every stress period/time step
  !<
  subroutine bnd_store_user_cond(this, nlist, rlist, condinput)
    ! -- modules
    use SimModule, only: store_error
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    integer(I4B), intent(in) :: nlist
    real(DP), dimension(:, :), pointer, contiguous, intent(in) :: rlist
    real(DP), dimension(:), pointer, contiguous, intent(inout) :: condinput
    ! -- local
    integer(I4B) :: l
    !
    ! -- store backup copy of conductance values
    do l = 1, nlist
      condinput(l) = rlist(2, l)
    end do
  end subroutine bnd_store_user_cond

  !> @ brief Read initial parameters for package
  !!
  !!  Read initial parameters for a boundary package. This method is not
  !!  needed for most boundary packages. The SFR package is an example of a
  !!  package that has overridden this method.
  !<
  subroutine bnd_read_initial_attr(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
  end subroutine bnd_read_initial_attr

  !> @ brief Read additional options for package
  !!
  !!  Read additional options for a boundary package. This method should
  !!  be overridden options in addition to the base options are implemented
  !!  in a boundary package.
  !<
  subroutine bnd_options(this, option, found)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    character(len=*), intent(inout) :: option !< option keyword string
    logical(LGP), intent(inout) :: found !< boolean indicating if the option was found
    !
    ! Return with found = .false.
    found = .false.
  end subroutine bnd_options

  !> @ brief Copy boundnames into boundnames_cst
  !!
  !!  boundnames_cst is an array of type(CharacterStringType),
  !!  which can be stored in the MemoryManager.
  !<
  subroutine copy_boundname(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    integer(I4B) :: i
    !
    ! copy from boundname to boundname_cst, which can be
    ! stored in the memory manager
    if (this%inamedbound /= 0) then
      do i = 1, size(this%boundname)
        this%boundname_cst(i) = this%boundname(i)
      end do
    end if
  end subroutine copy_boundname

  !> @ brief Setup output table for package
  !!
  !!  Setup output table for a boundary package that is used to output
  !!  package to model flow terms to the model listing file.
  !<
  subroutine pak_setup_outputtab(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    character(len=LINELENGTH) :: title
    character(len=LINELENGTH) :: text
    integer(I4B) :: ntabcol
    !
    ! -- allocate and initialize the output table
    if (this%iprflow /= 0) then
      !
      ! -- dimension table
      ntabcol = 3
      if (this%inamedbound > 0) then
        ntabcol = ntabcol + 1
      end if
      !
      ! -- initialize the output table object
      title = trim(adjustl(this%text))//' PACKAGE ('//trim(this%packName)// &
              ') FLOW RATES'
      call table_cr(this%outputtab, this%packName, title)
      call this%outputtab%table_df(this%maxbound, ntabcol, this%iout, &
                                   transient=.TRUE.)
      text = 'NUMBER'
      call this%outputtab%initialize_column(text, 10, alignment=TABCENTER)
      text = 'CELLID'
      call this%outputtab%initialize_column(text, 20, alignment=TABLEFT)
      text = 'RATE'
      call this%outputtab%initialize_column(text, 15, alignment=TABCENTER)
      if (this%inamedbound > 0) then
        text = 'NAME'
        call this%outputtab%initialize_column(text, LENBOUNDNAME, &
                                              alignment=TABLEFT)
      end if
    end if
  end subroutine pak_setup_outputtab

  !> @ brief Define the list label for the package
  !!
  !!  Method defined the list label for the boundary package. This method
  !!  needs to be overridden by each boundary package.
  !<
  subroutine define_listlabel(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
  end subroutine define_listlabel

  ! -- Procedures related to observations

  !> @brief Determine if observations are supported.
  !!
  !!  Function to determine if observations are supported by the boundary package.
  !!  By default, observations are not supported. This method should be overridden
  !!  if observations are supported in a boundary package.
  !!
  !!  @return  supported       boolean indicating if observations are supported
  !<
  function bnd_obs_supported(this) result(supported)
    ! -- return variable
    logical(LGP) :: supported !< boolean indicating if observations are supported
    ! -- dummy
    class(BndType) :: this !< BndType object
    !
    ! -- initialize return variables
    supported = .false.
  end function bnd_obs_supported

  !> @brief Define the observation types available in the package
  !!
  !!  Method to define the observation types available in a boundary
  !!  package. This method should be overridden if observations are
  !!  supported in a boundary package.
  !<
  subroutine bnd_df_obs(this)
    !
    ! -- dummy
    class(BndType) :: this !< BndType object
    !
    ! -- do nothing here. Override as needed.
  end subroutine bnd_df_obs

  !> @brief Read and prepare observations for a package
  !!
  !! Method to read and prepare observations for a boundary package
  !! This method should not need to be overridden for most boundary
  !! packages.
  !<
  subroutine bnd_rp_obs(this)
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: j
    class(ObserveType), pointer :: obsrv => null()
    character(len=LENBOUNDNAME) :: bname
    logical(LGP) :: jfound
    !
    if (.not. this%bnd_obs_supported()) return
    !
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      !
      ! -- indxbnds needs to be reset each stress period because
      !    list of boundaries can change each stress period.
      call obsrv%ResetObsIndex()
      obsrv%BndFound = .false.
      !
      bname = obsrv%FeatureName
      if (bname /= '') then
        !
        ! -- Observation location(s) is(are) based on a boundary name.
        !    Iterate through all boundaries to identify and store
        !    corresponding index(indices) in bound array.
        jfound = .false.
        do j = 1, this%nbound
          if (this%boundname(j) == bname) then
            jfound = .true.
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(j)
          end if
        end do
      else
        !
        ! -- Observation location is a single node number
        jfound = .false.
        jloop: do j = 1, this%nbound
          if (this%nodelist(j) == obsrv%NodeNumber) then
            jfound = .true.
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(j)
          end if
        end do jloop
      end if
    end do
    !
    if (count_errors() > 0) then
      call store_error_unit(this%inunit)
    end if
  end subroutine bnd_rp_obs

  !> @brief Save observations for the package
  !!
  !! Method to save simulated values for the boundary package.
  !! This method will need to be overridden for boundary packages
  !! with more observations than the calculate flow term (simvals)
  !! and to-mover.
  !<
  subroutine bnd_bd_obs(this)
    ! -- dummy
    class(BndType) :: this !< BndType object
    ! -- local
    integer(I4B) :: i
    integer(I4B) :: n
    real(DP) :: v
    type(ObserveType), pointer :: obsrv => null()
    !
    ! -- clear the observations
    call this%obs%obs_bd_clear()
    !
    ! -- Save simulated values for all of package's observations.
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      if (obsrv%BndFound) then
        do n = 1, obsrv%indxbnds_count
          if (obsrv%ObsTypeId == 'TO-MVR') then
            if (this%imover == 1) then
              v = this%pakmvrobj%get_qtomvr(obsrv%indxbnds(n))
              if (v > DZERO) then
                v = -v
              end if
            else
              v = DNODATA
            end if
          else
            v = this%simvals(obsrv%indxbnds(n))
          end if
          call this%obs%SaveOneSimval(obsrv, v)
        end do
      else
        call this%obs%SaveOneSimval(obsrv, DNODATA)
      end if
    end do
  end subroutine bnd_bd_obs

  !> @brief Output observations for the package
  !!
  !! Method to output simulated values for the boundary package.
  !! This method should not need to be overridden.
  !<
  subroutine bnd_ot_obs(this)
    ! -- dummy
    class(BndType) :: this !< BndType object
    !
    ! -- call the observation output method
    call this%obs%obs_ot()
  end subroutine bnd_ot_obs

  ! -- Procedures related to time series

  !> @brief Assign time series links for the package
  !!
  !! Assign the time series links for the boundary package. This
  !! method will need to be overridden for boundary packages that
  !! support time series.
  !<
  subroutine bnd_rp_ts(this)
    ! -- dummy
    class(BndType), intent(inout) :: this
  end subroutine bnd_rp_ts

  ! -- Procedures related to casting

  !> @brief Cast as a boundary type
  !!
  !!  Subroutine to cast an object as a boundary package type.
  !<
  function CastAsBndClass(obj) result(res)
    class(*), pointer, intent(inout) :: obj !< input object
    class(BndType), pointer :: res !< output class of type BndType
    !
    ! -- initialize res
    res => null()
    !
    ! -- make sure obj is associated
    if (.not. associated(obj)) return
    !
    ! -- point res to obj
    select type (obj)
    class is (BndType)
      res => obj
    end select
  end function CastAsBndClass

  !> @brief Add boundary to package list
  !!
  !!  Subroutine to add a boundary package to a package list.
  !<
  subroutine AddBndToList(list, bnd)
    ! -- dummy
    type(ListType), intent(inout) :: list !< package list
    class(BndType), pointer, intent(inout) :: bnd !< boundary package
    ! -- local
    class(*), pointer :: obj
    !
    obj => bnd
    call list%Add(obj)
  end subroutine AddBndToList

  !> @brief Get boundary from package list
  !!
  !!  Function to get a boundary package from a package list.
  !!
  !!  @return  res  boundary package object
  !<
  function GetBndFromList(list, idx) result(res)
    ! -- dummy
    type(ListType), intent(inout) :: list !< package list
    integer(I4B), intent(in) :: idx !< package number
    class(BndType), pointer :: res !< boundary package idx
    ! -- local
    class(*), pointer :: obj
    !
    ! -- get the package from the list
    obj => list%GetItem(idx)
    res => CastAsBndClass(obj)
  end function GetBndFromList

  !> @brief Save and/or print flows for a package
  !!
  !!  Subroutine to save and/or print package flows to a model to a
  !!  binary cell-by-cell flow file and the model listing file.
  !<
  subroutine save_print_model_flows(icbcfl, ibudfl, icbcun, iprflow, &
                                    outputtab, nbound, nodelist, flow, ibound, &
                                    title, text, ipakcb, dis, naux, textmodel, &
                                    textpackage, dstmodel, dstpackage, &
                                    auxname, auxvar, iout, inamedbound, &
                                    boundname, imap)
    ! -- modules
    use TdisModule, only: kstp, kper
    ! -- dummy
    integer(I4B), intent(in) :: icbcfl !< flag indicating if the flow should be saved to the binary cell-by-cell flow file
    integer(I4B), intent(in) :: ibudfl !< flag indicating if the flow should be saved or printed
    integer(I4B), intent(in) :: icbcun !< file unit number for the binary cell-by-cell file
    integer(I4B), intent(in) :: iprflow !< print flows to list file
    type(TableType), pointer, intent(inout) :: outputtab !< output table object
    integer(I4B), intent(in) :: nbound !< number of boundaries this stress period
    integer(I4B), dimension(:), contiguous, intent(in) :: nodelist !< boundary node list
    real(DP), dimension(:), contiguous, intent(in) :: flow !< boundary flow terms
    integer(I4B), dimension(:), contiguous, intent(in) :: ibound !< ibound array for the model
    character(len=*), intent(in) :: title !< title for the output table
    character(len=*), intent(in) :: text !< flow term description
    integer(I4B), intent(in) :: ipakcb !< flag indicating if flows will be saved
    class(DisBaseType), intent(in) :: dis !< model discretization object
    integer(I4B), intent(in) :: naux !< number of aux variables
    character(len=*), intent(in) :: textmodel !< model name
    character(len=*), intent(in) :: textpackage !< package name
    character(len=*), intent(in) :: dstmodel !< mover destination model
    character(len=*), intent(in) :: dstpackage !< mover destination package
    character(len=*), dimension(:), intent(in) :: auxname !< aux variable name
    real(DP), dimension(:, :), intent(in) :: auxvar !< aux variable
    integer(I4B), intent(in) :: iout !< model listing file unit
    integer(I4B), intent(in) :: inamedbound !< flag indicating if boundnames are defined for the boundary entries
    character(len=LENBOUNDNAME), dimension(:), contiguous :: boundname !< bound names
    integer(I4B), dimension(:), optional, intent(in) :: imap !< mapping array
    ! -- local
    character(len=20) :: nodestr
    integer(I4B) :: nodeu
    integer(I4B) :: maxrows
    integer(I4B) :: i
    integer(I4B) :: node
    integer(I4B) :: n2
    integer(I4B) :: ibinun
    integer(I4B) :: nboundcount
    real(DP) :: rrate
    real(DP), dimension(naux) :: auxrow
    ! -- for observations
    character(len=LENBOUNDNAME) :: bname
    !
    ! -- set table kstp and kper
    if (iprflow /= 0) then
      call outputtab%set_kstpkper(kstp, kper)
    end if
    !
    ! -- set maxrows
    maxrows = 0
    if (ibudfl /= 0 .and. iprflow /= 0) then
      do i = 1, nbound
        node = nodelist(i)
        if (node > 0) then
          maxrows = maxrows + 1
        end if
      end do
      if (maxrows > 0) then
        call outputtab%set_maxbound(maxrows)
      end if
      call outputtab%set_title(title)
    end if
    !
    ! -- Set unit number for binary output
    if (ipakcb < 0) then
      ibinun = icbcun
    else if (ipakcb == 0) then
      ibinun = 0
    else
      ibinun = ipakcb
    end if
    if (icbcfl == 0) then
      ibinun = 0
    end if
    !
    ! -- If cell-by-cell flows will be saved as a list, write header.
    if (ibinun /= 0) then
      !
      ! -- Count nbound as the number of entries with node > 0
      !    SFR, for example, can have a 'none' connection, which
      !    means it should be excluded from budget file
      nboundcount = 0
      do i = 1, nbound
        node = nodelist(i)
        if (node > 0) nboundcount = nboundcount + 1
      end do
      call dis%record_srcdst_list_header(text, textmodel, textpackage, &
                                         dstmodel, dstpackage, naux, &
                                         auxname, ibinun, nboundcount, iout)
    end if
    !
    ! -- If no boundaries, skip flow calculations.
    if (nbound > 0) then
      !
      ! -- Loop through each boundary calculating flow.
      do i = 1, nbound
        node = nodelist(i)
        ! -- assign boundary name
        if (inamedbound > 0) then
          bname = boundname(i)
        else
          bname = ''
        end if
        !
        ! -- If cell is no-flow or constant-head, then ignore it.
        rrate = DZERO
        if (node > 0) then
          !
          ! -- Use simval, which was calculated in cq()
          rrate = flow(i)
          !
          ! -- Print the individual rates if the budget is being printed
          !    and PRINT_FLOWS was specified (iprflow < 0).  Rates are
          !    printed even if ibound < 1.
          if (ibudfl /= 0) then
            if (iprflow /= 0) then
              !
              ! -- set nodestr and write outputtab table
              nodeu = dis%get_nodeuser(node)
              call dis%nodeu_to_string(nodeu, nodestr)
              call outputtab%print_list_entry(i, trim(adjustl(nodestr)), &
                                              rrate, bname)
            end if
          end if
          !
          ! -- If saving cell-by-cell flows in list, write flow
          if (ibinun /= 0) then
            n2 = i
            if (present(imap)) n2 = imap(i)
            if (naux > 0) then
              auxrow(:) = auxvar(:, i)
            end if
            call dis%record_mf6_list_entry(ibinun, node, n2, rrate, naux, &
                                           auxrow, olconv2=.FALSE.)
          end if
        end if
        !
      end do
      if (ibudfl /= 0) then
        if (iprflow /= 0) then
          write (iout, '(1x)')
        end if
      end if

    end if
  end subroutine save_print_model_flows

  !> @brief Activate viscosity terms
  !!
  !! Method to activate addition of viscosity terms when package type
  !! is DRN, GHB, or RIV (method not needed by other packages at this point)
  !<
  subroutine bnd_activate_viscosity(this)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate
    ! -- dummy
    class(BndType), intent(inout) :: this !< BndType object
    ! -- local
    integer(I4B) :: i
    !
    ! -- Set ivsc and reallocate viscratios to be of size MAXBOUND
    this%ivsc = 1
    !
    ! -- Allocate array for storing user-specified conductances
    !    modified by updated viscosity values
    call mem_reallocate(this%condinput, this%maxbound, 'CONDINPUT', &
                        this%memoryPath)
    do i = 1, this%maxbound
      this%condinput(i) = DZERO
    end do
    !
    ! -- Notify user via listing file viscosity accounted for by standard
    !    boundary package.
    write (this%iout, '(/1x,a,a)') 'VISCOSITY ACTIVE IN ', &
      trim(this%filtyp)//' PACKAGE CALCULATIONS: '//trim(adjustl(this%packName))
  end subroutine bnd_activate_viscosity

end module BndModule
