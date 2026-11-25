!> @brief This module contains the CSUB package methods
!!
!! This module contains the methods used to add the effects of elastic
!! skeletal storage, compaction, and subsidence on the groundwater flow
!! equation. The contribution of elastic skelatal, inelastic and elastic
!! interbed storage and water compressibility can be represented.
!!
!! This file holds the type, the interfaces, the create, allocate and read,
!! read and prepare, advance, and fill routines, and the shared helpers. The
!! other procedures are in submodules/gwf-csub-*.f90:
!!
!! - memory: allocation and deallocation
!! - input: options, dimensions, griddata, and packagedata
!! - stress: initial states and coarse-grained stresses
!! - coarse: coarse-grained storage
!! - nodelay: no-delay interbed storage
!! - delay: delay interbeds
!! - budget: convergence check, flows, and budget
!! - output: final processing and dependent variable output
!! - obs: observations
!!
!<
module GwfCsubModule
  use KindModule, only: I4B, DP, LGP
  use ConstantsModule, only: DPREC, DZERO, DEM20, DEM15, DEM10, DEM8, DEM7, &
                             DEM6, DEM5, DEM4, DEM3, DP9, DHALF, DEM1, DONE, &
                             DTWO, DTHREE, DGRAVITY, DTEN, DHUNDRED, DNODATA, &
                             DHNOFLO, LENFTYPE, LENPACKAGENAME, LENMEMPATH, &
                             LINELENGTH, LENBOUNDNAME, NAMEDBOUNDFLAG, &
                             LENBUDTXT, LENAUXNAME, LENPAKLOC, LENLISTLABEL, &
                             TABLEFT, TABCENTER, TABRIGHT, &
                             TABSTRING, TABUCSTRING, TABINTEGER, TABREAL
  use MemoryHelperModule, only: create_mem_path
  use MathUtilModule, only: is_close
  use MessageModule, only: write_message
  use SmoothingModule, only: sQuadraticSaturation, &
                             sQuadraticSaturationDerivative, &
                             sQuadratic0sp, &
                             sQuadratic0spDerivative
  use NumericalPackageModule, only: NumericalPackageType
  use ObserveModule, only: ObserveType
  use ObsModule, only: ObsType, obs_cr
  use GeomUtilModule, only: get_node
  use InputOutputModule, only: extract_idnum_or_bndname
  use BaseDisModule, only: DisBaseType
  use SimModule, only: count_errors, store_error, store_warning, &
                       store_error_filename
  use SimVariablesModule, only: errmsg, warnmsg
  use SortModule, only: selectn
  !
  use TableModule, only: TableType, table_cr
  use BudgetModule, only: BudgetType
  use CompactionUtilModule, only: porosity_to_void_ratio, compact_layer, &
                                  adjust_eff_stress, storage_factor, &
                                  storage_newton_term, wcomp_newton_term
  !
  use MatrixBaseModule
  !
  implicit none
  !
  private
  public :: csub_cr
  public :: GwfCsubType
  !
  ! -- local parameter
  real(DP), parameter :: dlog10es = 0.4342942_DP !< derivative of the log of effective stress
  real(DP), parameter :: stressfloor = DEM3 !< effective-stress regularization floor (fraction of geostatic stress)
  !
  ! CSUB type
  type, extends(NumericalPackageType) :: GwfCsubType
    ! -- characters scalars
    character(len=LENLISTLABEL), pointer :: listlabel => null() !< title of table written for RP
    character(len=LENMEMPATH), pointer :: stoMemPath => null() !< memory path of storage package
    ! -- character arrays
    character(len=LENBOUNDNAME), dimension(:), &
      pointer, contiguous :: boundname => null() !< vector of boundnames
    character(len=LENAUXNAME), dimension(:), &
      pointer, contiguous :: auxname => null() !< vector of auxname
    ! -- logical scalars
    logical(LGP), pointer :: lhead_based => null() !< logical variable indicating if head-based solution
    ! -- integer scalars
    integer(I4B), pointer :: istounit => null() !< unit number of storage package
    integer(I4B), pointer :: istrict_stress => null() !< 1 terminates on negative effective stress, 0 regularizes (default)
    integer(I4B), pointer :: nreg_ts => null() !< count of time steps with effective-stress regularization
    integer(I4B), pointer :: istrainib => null() !< unit number of interbed strain output
    integer(I4B), pointer :: istrainsk => null() !< unit number of coarse-grained strain output
    integer(I4B), pointer :: ioutcomp => null() !< unit number for cell-by-cell compaction output
    integer(I4B), pointer :: ioutcompi => null() !< unit number for cell-by-cell inelastic compaction output
    integer(I4B), pointer :: ioutcompe => null() !< unit number for cell-by-cell elastic compaction output
    integer(I4B), pointer :: ioutcompib => null() !< unit number for cell-by-cell interbed compaction output
    integer(I4B), pointer :: ioutcomps => null() !< unit number for cell-by-cell coarse-grained compaction output
    integer(I4B), pointer :: ioutzdisp => null() !< unit number for z-displacement output
    integer(I4B), pointer :: ipakcsv => null() !< unit number for csv output
    integer(I4B), pointer :: iupdatematprop => null() !< flag indicating if material properties will be updated
    integer(I4B), pointer :: istoragec => null() !< flag indicating specific storage coefficient will be specified
    integer(I4B), pointer :: icellf => null() !< flag indicating cell fractions will be specified
    integer(I4B), pointer :: ispecified_pcs => null() !< flag indicating preconsolidation state is specified (not relative)
    integer(I4B), pointer :: ispecified_dbh => null() !< flag indicating delay bed head is specified (not relative)
    integer(I4B), pointer :: inamedbound => null() !< flag to read boundnames
    integer(I4B), pointer :: iconvchk => null() !< flag indicating if a final convergence check will be made
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    integer(I4B), pointer :: ninterbeds => null() !< number of interbeds
    integer(I4B), pointer :: maxsig0 => null() !< maximum number of cells with specified sig0 values
    integer(I4B), pointer :: nbound => null() !< number of boundaries for current stress period
    integer(I4B), pointer :: iscloc => null() !< bound column to scale with SFAC
    integer(I4B), pointer :: iauxmultcol => null() !< column to use as multiplier for column iscloc
    integer(I4B), pointer :: ndelaycells => null() !< number of cells in delay interbeds
    integer(I4B), pointer :: ndelaybeds => null() !< number of delay interbeds
    integer(I4B), pointer :: initialized => null() !< flag indicating if the initial stresses have been initialized
    integer(I4B), pointer :: ieslag => null() !< flag indicating if the effective stress is lagged
    integer(I4B), pointer :: ipch => null() !< flag indicating if initial precosolidation value is a head
    integer(I4B), pointer :: iupdatestress => null() !< flag indicating if the geostatic stress is active
    ! -- real scalars
    real(DP), pointer :: epsilon => null() !< epsilon for stress smoothing
    real(DP), pointer :: cc_crit => null() !< convergence criteria for csub-gwf convergence check
    real(DP), pointer :: gammaw => null() !< product of fluid density, and gravity
    real(DP), pointer :: beta => null() !< water compressibility
    real(DP), pointer :: brg => null() !< product of gammaw and water compressibility
    real(DP), pointer :: satomega => null() !< newton-raphson saturation omega
    real(DP), pointer :: pcsomega => null() !< elastic<->inelastic switch smoothing window (fraction of pcs; 0 = hard switch)
    ! -- integer pointer to storage package variables
    integer(I4B), pointer :: gwfiss => NULL() !< pointer to model iss flag
    integer(I4B), pointer :: gwfiss0 => NULL() !< iss flag for last stress period
    ! -- integer arrays
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    integer(I4B), dimension(:), pointer, contiguous :: stoiconv => null() !< pointer to iconvert in storage
    ! -- real arrays
    real(DP), dimension(:), pointer, contiguous :: stoss => null() !< pointer to ss in storage
    real(DP), dimension(:), pointer, contiguous :: buff => null() !< buff array
    real(DP), dimension(:), pointer, contiguous :: buffusr => null() !< buffusr array
    integer, dimension(:), pointer, contiguous :: nodelist => null() !< reduced node that the interbed is attached to
    integer, dimension(:), pointer, contiguous :: unodelist => null() !< user node that the interbed is attached to
    !
    ! -- coarse-grained storage variables
    real(DP), dimension(:), pointer, contiguous :: sgm => null() !< specific gravity moist sediments
    real(DP), dimension(:), pointer, contiguous :: sgs => null() !< specific gravity saturated sediments
    real(DP), dimension(:), pointer, contiguous :: cg_ske_cr => null() !< coarse-grained specified storage
    real(DP), dimension(:), pointer, contiguous :: cg_gs => null() !< geostatic stress for a cell
    real(DP), dimension(:), pointer, contiguous :: cg_es => null() !< coarse-grained (aquifer) effective stress
    real(DP), dimension(:), pointer, contiguous :: cg_es0 => null() !< coarse-grained (aquifer) effective stress for the previous time step
    real(DP), dimension(:), pointer, contiguous :: cg_pcs => null() !< coarse-grained (aquifer) preconsolidation stress
    real(DP), dimension(:), pointer, contiguous :: cg_comp => null() !< coarse-grained (aquifer) incremental compaction
    real(DP), dimension(:), pointer, contiguous :: cg_tcomp => null() !< coarse-grained (aquifer) total compaction
    real(DP), dimension(:), pointer, contiguous :: cg_stor => null() !< coarse-grained (aquifer) storage
    real(DP), dimension(:), pointer, contiguous :: cg_ske => null() !< coarse-grained (aquifer) elastic storage coefficient
    real(DP), dimension(:), pointer, contiguous :: cg_sk => null() !< coarse-grained (aquifer) first storage coefficient
    real(DP), dimension(:), pointer, contiguous :: cg_thickini => null() !< initial coarse-grained (aquifer) thickness
    real(DP), dimension(:), pointer, contiguous :: cg_thetaini => null() !< initial coarse-grained (aquifer) porosity
    real(DP), dimension(:), pointer, contiguous :: cg_thick => null() !< current coarse-grained (aquifer) thickness
    real(DP), dimension(:), pointer, contiguous :: cg_thick0 => null() !< previous coarse-grained (aquifer) thickness
    real(DP), dimension(:), pointer, contiguous :: cg_theta => null() !< current coarse-grained (aquifer) porosity
    real(DP), dimension(:), pointer, contiguous :: cg_theta0 => null() !< previous coarse-grained (aquifer) porosity
    !
    ! -- cell storage variables
    real(DP), dimension(:), pointer, contiguous :: cell_wcstor => null() !< cell water compressibility storage
    real(DP), dimension(:), pointer, contiguous :: cell_thick => null() !< cell compressible material thickness
    !
    ! -- interbed variables
    integer(I4B), dimension(:), pointer, contiguous :: idelay => null() !< delay interbed flag - 0 = nodelay, > 0 = delay
    integer(I4B), dimension(:), pointer, contiguous :: ielastic => null() !< elastic interbed equation - 0 = inelastic and elastic, > 0 = elastic
    integer(I4B), dimension(:), pointer, contiguous :: iconvert => null() !< convertible cell flag - 0 = elastic, > 0 = inelastic
    integer(I4B), dimension(:), pointer, contiguous :: iconvert_nr => null() !< iconvert at the last Newton-Raphson term (-1 = not yet set)
    real(DP), dimension(:), pointer, contiguous :: ci => null() !< compression index
    real(DP), dimension(:), pointer, contiguous :: rci => null() !< recompression index
    real(DP), dimension(:), pointer, contiguous :: pcs => null() !< preconsolidation stress
    real(DP), dimension(:), pointer, contiguous :: rnb => null() !< interbed system material factor
    real(DP), dimension(:), pointer, contiguous :: kv => null() !< vertical hydraulic conductivity of interbed
    real(DP), dimension(:), pointer, contiguous :: h0 => null() !< initial head in interbed
    real(DP), dimension(:), pointer, contiguous :: comp => null() !< interbed incremental compaction
    real(DP), dimension(:), pointer, contiguous :: tcomp => null() !< total interbed compaction
    real(DP), dimension(:), pointer, contiguous :: tcompi => null() !< total inelastic interbed compaction
    real(DP), dimension(:), pointer, contiguous :: tcompe => null() !< total elastic interbed compaction
    real(DP), dimension(:), pointer, contiguous :: storagee => null() !< elastic storage
    real(DP), dimension(:), pointer, contiguous :: storagei => null() !< inelastic storage
    real(DP), dimension(:), pointer, contiguous :: ske => null() !< elastic storage coefficient
    real(DP), dimension(:), pointer, contiguous :: sk => null() !< first storage coefficient
    real(DP), dimension(:), pointer, contiguous :: thickini => null() !< initial interbed thickness
    real(DP), dimension(:), pointer, contiguous :: thetaini => null() !< initial interbed theta
    real(DP), dimension(:), pointer, contiguous :: thick => null() !< current interbed thickness
    real(DP), dimension(:), pointer, contiguous :: thick0 => null() !< previous interbed thickness
    real(DP), dimension(:), pointer, contiguous :: theta => null() !< current interbed porosity
    real(DP), dimension(:), pointer, contiguous :: theta0 => null() !< previous interbed porosity
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< auxiliary variable array
    !
    ! -- delay interbed
    integer(I4B), dimension(:), pointer, contiguous :: idb_nconv_count => null() !< non-convertible count of interbeds with heads below delay cell top
    integer(I4B), dimension(:, :), pointer, contiguous :: idbconvert => null() !0 = elastic, > 0 = inelastic
    real(DP), dimension(:), pointer, contiguous :: dbdhmax => null() !< delay bed maximum head change
    real(DP), dimension(:, :), pointer, contiguous :: dbz => null() !< delay bed cell z
    real(DP), dimension(:, :), pointer, contiguous :: dbrelz => null() !< delay bed cell z relative to znode
    real(DP), dimension(:, :), pointer, contiguous :: dbh => null() !< delay bed cell h
    real(DP), dimension(:, :), pointer, contiguous :: dbh0 => null() !< delay bed cell previous h
    real(DP), dimension(:, :), pointer, contiguous :: dbgeo => null() !< delay bed cell geostatic stress
    real(DP), dimension(:, :), pointer, contiguous :: dbes => null() !< delay bed cell effective stress
    real(DP), dimension(:, :), pointer, contiguous :: dbes0 => null() !< delay bed cell previous effective stress
    real(DP), dimension(:, :), pointer, contiguous :: dbpcs => null() !< delay bed cell preconsolidation stress
    real(DP), dimension(:), pointer, contiguous :: dbflowtop => null() !< delay bed flow through interbed top
    real(DP), dimension(:), pointer, contiguous :: dbflowbot => null() !< delay bed flow through interbed bottom
    real(DP), dimension(:, :), pointer, contiguous :: dbdzini => null() !< initial delay bed cell thickness
    real(DP), dimension(:, :), pointer, contiguous :: dbthetaini => null() !< initial delay bed cell porosity
    real(DP), dimension(:, :), pointer, contiguous :: dbdz => null() !< delay bed dz
    real(DP), dimension(:, :), pointer, contiguous :: dbdz0 => null() !< delay bed previous dz
    real(DP), dimension(:, :), pointer, contiguous :: dbtheta => null() !< delay bed cell porosity
    real(DP), dimension(:, :), pointer, contiguous :: dbtheta0 => null() !< delay bed cell previous porosity
    real(DP), dimension(:, :), pointer, contiguous :: dbcomp => null() !< delay bed incremental compaction
    real(DP), dimension(:, :), pointer, contiguous :: dbtcomp => null() !< delay bed total interbed compaction
    !
    ! -- delay interbed solution arrays
    real(DP), dimension(:), pointer, contiguous :: dbal => null() !< delay bed lower diagonal
    real(DP), dimension(:), pointer, contiguous :: dbad => null() !< delay bed diagonal
    real(DP), dimension(:), pointer, contiguous :: dbau => null() !< delay bed upper diagonal
    real(DP), dimension(:), pointer, contiguous :: dbrhs => null() !< delay bed right hand side
    real(DP), dimension(:), pointer, contiguous :: dbdh => null() !< delay bed dh
    real(DP), dimension(:), pointer, contiguous :: dbaw => null() !< delay bed work vector
    !
    ! -- period data
    integer(I4B), dimension(:), pointer, contiguous :: nodelistsig0 => null() !< vector of reduced node numbers
    real(DP), dimension(:), pointer, contiguous :: sig0 => null() !< array of package specific boundary numbers
    !
    ! -- observation data
    integer(I4B), pointer :: inobspkg => null() !< unit number for obs package
    type(ObsType), pointer :: obs => null() !< observation package
    !
    ! -- table objects
    type(TableType), pointer :: inputtab => null() !< table for input variables
    type(TableType), pointer :: outputtab => null() !< table for output variables
    type(TableType), pointer :: pakcsvtab => null() !< table for csv output

  contains
    !> @brief Define the list label for the package
    !!
    !!  Method defined the list label for the CSUB package. The list label is
    !!  the heading that is written to iout when PRINT_INPUT option is used.
    !!
    !<
    procedure :: define_listlabel
    !> @brief Source options for package
    !!
    !!  Source options for CSUB package.
    !!
    !<
    procedure :: source_options
    !> @brief log options for package
    !!
    !!  log options block for CSUB package.
    !!
    !<
    procedure :: log_options
    procedure :: csub_ar
    !> @brief Deallocate package memory
    !!
    !!  Deallocate CSUB package scalars and arrays.
    !!
    !<
    procedure :: csub_da
    procedure :: csub_rp
    procedure :: csub_ad
    procedure :: csub_fc
    procedure :: csub_fn
    !> @brief Final convergence check
    !!
    !! Final convergence check for the CSUB package. The final convergence
    !! check is only required when the simulation includes delay interbeds.
    !! The final convergence check compares the sum of water contributed
    !! by storage and water compressibility in the delay bed to the fluid
    !! exchange between the delay interbed and the gwf cell.
    !<
    procedure :: csub_cc
    !> @brief Calculate flows for package
    !!
    !!  Flow calculation for the CSUB package components. Components include
    !!  coarse-grained storage, delay and no-delay interbeds, and water
    !!  compressibility.
    !!
    !<
    procedure :: csub_cq
    !> @brief Model budget calculation for package
    !!
    !!  Budget calculation for the CSUB package components. Components include
    !!  coarse-grained storage, delay and no-delay interbeds, and water
    !!  compressibility.
    !<
    procedure :: csub_bd
    !> @brief Save model flows for package
  !!
  !!  Save cell-by-cell budget terms for the CSUB package.
  !!
    !<
    procedure :: csub_save_model_flows
    !> @brief Save and print dependent values for package
  !!
  !!  Method saves cell-by-cell compaction and z-displacement terms. The method
  !!  also calls the method to process observation output.
  !!
    !<
    procedure :: csub_ot_dv
    !> @brief Final processing for package
    !!
    !!  Final processing for the  CSUB package. This method generates the final
    !!  strain tables that are output so that the user can evaluate if calculated
    !!  strain rates in coarse-grained sediments and interbeds exceed 1 percent.
    !!
    !<
    procedure :: csub_fp
    !> @brief Source dimensions for package
    !!
    !!  Read the number of interbeds and maximum number of cells with a specified
    !!  overlying geostatic stress.
    !!
    !<
    procedure :: source_dimensions => csub_source_dimensions
    !> @brief Allocate scalars
    !!
    !! Allocate and initialize scalars for the CSUB package. The base model
    !! allocate scalars method is also called.
    !!
    !<
    procedure, private :: csub_allocate_scalars
    !> @brief Allocate package arrays
    !!
    !!  Allocate and initialize CSUB package arrays.
    !!
    !<
    procedure, private :: csub_allocate_arrays
    !> @brief Source griddata for package
    !<
    procedure, private :: csub_source_griddata
    !> @brief source packagedata for package
    !!
    !!  Read delay and no-delay interbed input data for the CSUB package. Method
    !!  also validates interbed input data.
    !!
    !<
    procedure, private :: csub_source_packagedata
    !> @brief Print packagedata
    !<
    procedure, private :: csub_print_packagedata
    !
    ! -- helper methods
    procedure, private :: csub_calc_znode
    procedure, private :: csub_calc_sat
    procedure, private :: csub_calc_sat_derivative
    procedure, private :: csub_calc_sfacts
    procedure, private :: csub_calc_interbed_thickness
    !> @brief Calculate the flow from delay interbed top or bottom
    !!
    !! Function to calculate the flow from across the top or bottom of
    !! a delay interbed.
    !!
    !! @return  q  flow across the top or bottom of a delay interbed
    !<
    procedure, private :: csub_calc_delay_flow
    !
    ! -- stress methods
    !> @brief Calculate the stress for model cells
    !!
    !!  Method calculates the geostatic stress, pressure head, and effective
    !!  stress at the bottom of each cell. The method also applies the overlying
    !!  geostatic stress (sig0) not represented in the model.
    !!
    !<
    procedure, private :: csub_cg_calc_stress
    !> @brief Check effective stress values
    !!
    !!  Method checks calculated effective stress values to ensure that
    !!  effective stress values are positive. An error condition and message
    !!  are issued if calculated effective stress values are less than a small
    !!  positive value (DEM6).
    !!
    !<
    procedure, private :: csub_cg_chk_stress
    !
    ! -- initial states
    !> @brief Set initial states for the package
    !!
    !!  Method sets the initial states for coarse-grained materials and fine-
    !!  grained sediments in the interbeds.
    !!
    !<
    procedure, private :: csub_set_initial_state
    !
    ! -- coarse-grained coarse-grained methods
    !> @brief Update coarse-grained material properties
    !!
    !!  Method updates coarse-grained material properties in a cell.
    !!
    !<
    procedure, private :: csub_cg_update
    !> @brief Calculate coarse-grained compaction in a cell
    !!
    !!  Method calculates coarse-grained compaction in a cell.
    !<
    procedure, private :: csub_cg_calc_comp
    !> @brief Calculate Sske for a cell
    !!
    !!  Method calculates Sske for coarse-grained materials in a cell.
    !<
    procedure, private :: csub_cg_calc_sske
    !> @brief Formulate the coefficients for coarse-grained materials
    !!
    !!  Method formulates the coefficient matrix and right-hand side terms
    !!  for coarse grained materials in a cell.
    !<
    procedure, private :: csub_cg_fc
    !> @brief Formulate coarse-grained Newton-Raphson terms
    !!
    !!  Method formulates the coefficient matrix and right-hand side terms
    !!  for coarse grained materials in a cell when using the Newton-Raphson
    !!  formulation.
    !<
    procedure, private :: csub_cg_fn
    !> @brief Formulate coarse-grained water compressibility coefficients
    !!
    !!  Method formulates the standard formulation coefficient matrix and
    !!  right-hand side terms for water compressibility in coarse-grained
    !!  sediments.
    !<
    procedure, private :: csub_cg_wcomp_fc
    !> @brief Formulate coarse-grained water compressibility coefficients
    !!
    !!  Method formulates the Newton-Raphson formulation coefficient matrix and
    !!  right-hand side terms for water compressibility in coarse-grained
    !!  sediments.
    !<
    procedure, private :: csub_cg_wcomp_fn
    !
    ! -- interbed methods
    !> @brief Formulate the coefficients for a interbed
    !!
    !!  Method formulates the coefficient matrix and right-hand side terms
    !!  for a interbed in a cell.
    !<
    procedure, private :: csub_interbed_fc
    !> @brief Formulate the coefficients for a interbed
    !!
    !!  Method formulates the Newton-Raphson formulation coefficient matrix and
    !!  right-hand side terms for a interbed in a cell.
    !<
    procedure, private :: csub_interbed_fn
    !
    ! -- no-delay interbed methods
    !> @brief Update no-delay material properties
    !!
    !!  Method updates no-delay material properties based on the current
    !!  compaction value.
    !!
    !<
    procedure, private :: csub_nodelay_update
    !> @brief Calculate no-delay interbed storage coefficients
    !!
    !!  Method calculates the skeletal storage coefficients for a no-delay
    !!  interbed. The method also calculates the contribution of the
    !!  no-delay interbed to the right-hand side of the groundwater flow
    !!  equation for the cell.
    !<
    procedure, private :: csub_nodelay_fc
    !> @brief Formulate no-delay interbed water compressibility coefficients
    !!
    !!  Method formulates the standard formulation coefficient matrix and
    !!  right-hand side terms for water compressibility in no-delay
    !!  interbeds.
    !<
    procedure, private :: csub_nodelay_wcomp_fc
    !> @brief Formulate no-delay interbed water compressibility coefficients
    !!
    !!  Method formulates the Newton-Raphson formulation coefficient matrix and
    !!  right-hand side terms for water compressibility in no-delay
    !!  interbeds.
    !<
    procedure, private :: csub_nodelay_wcomp_fn
    !> @brief Calculate no-delay interbed compaction
    !!
    !!  Method calculates the compaction for a no-delay interbed. The method
    !!  also calculates the storage coefficients for the no-delay interbed.
    !<
    procedure, private :: csub_nodelay_calc_comp
    !
    ! -- delay interbed methods
    !> @brief Calculate delay interbed saturation
    !!
    !! Method to calculate the saturation in a delay interbed cell.
    !<
    procedure, private :: csub_delay_calc_sat
    !> @brief Calculate the delay interbed cell saturation derivative
    !!
    !! Function to calculate the derivative of the saturation with
    !! respect to the current head in delay interbed cell n.
    !!
    !! @return      satderv              derivative of saturation
    !<
    procedure, private :: csub_delay_calc_sat_derivative
    !> @brief Calculate delay interbed znode and z relative to interbed center
    !!
    !! Method to calculate the initial center of each delay interbed cell,
    !! assuming the delay bed head is equal to the top of the delay interbed.
    !! The method also calculates the distance of the center of each delay
    !! bed cell from the center of the delay interbed (z_offset) that is used
    !! to calculate average skeletal specific storage values for a delay interbed
    !! centered on the center of the saturated thickness for a cell.
    !!
    !<
    procedure, private :: csub_delay_init_zcell
    !> @brief Calculate delay interbed stress values
    !!
    !! Method to calculate the geostatic and effective stress in delay interbed
    !! cells using the passed the current head value in a cell.
    !!
    !<
    procedure, private :: csub_delay_calc_stress
    !> @brief Calculate delay interbed cell storage coefficients
    !<
    procedure, private :: csub_delay_calc_ssksske
    !> @brief Calculate delay interbed compaction
    !!
    !! Method to calculate the compaction in a delay interbed.
    !<
    procedure, private :: csub_delay_calc_comp
    !> @brief Update delay interbed material properties
    !!
    !! Method to update the thickness and porosity of each delay interbed cell.
    !!
    !<
    procedure, private :: csub_delay_update
    !> @brief Calculate delay interbed storage change
    !!
    !! Method to calculate the storage change in a delay interbed.
    !<
    procedure, private :: csub_delay_calc_dstor
    !> @brief Calculate delay interbed water compressibility
    !!
    !! Method to calculate the change in water compressibility in a delay
    !! interbed.
    !<
    procedure, private :: csub_delay_calc_wcomp
    !> @brief Calculate delay interbed contribution to the cell
    !!
    !! Method to calculate the coefficients to calculate the delay interbed
    !! contribution to a cell. The product of hcof* h - rhs equals the
    !! delay contribution to the cell
    !<
    procedure, private :: csub_delay_fc
    !> @brief Solve delay interbed continuity equation
    !!
    !! Method to calculate solve the delay interbed continuity equation for a
    !! delay interbed. The method encapsulates the non-linear loop and calls the
    !! linear solution.
    !!
    !<
    procedure, private :: csub_delay_sln
    !> @brief Assemble delay interbed coefficients
    !!
    !! Method to assemble matrix and right-hand side coefficients for a delay
    !! interbed. The method calls the appropriate standard or Newton-Raphson
    !! assembly routines and fills all of the entries for a delay interbed.
    !!
    !<
    procedure, private :: csub_delay_assemble
    !> @brief Assemble delay interbed standard formulation coefficients
    !!
    !! Method to assemble standard formulation matrix and right-hand side
    !! coefficients for a delay interbed.
    !!
    !<
    procedure, private :: csub_delay_assemble_fc
    !> @brief Assemble delay interbed Newton-Raphson formulation coefficients
    !!
    !! Method to assemble Newton-Raphson formulation matrix and right-hand side
    !! coefficients for a delay interbed.
    !!
    !<
    procedure, private :: csub_delay_assemble_fn
    !> @brief Check delay interbed head
    !!
    !! Method to determine if the delay interbed head in any delay cell
    !! in a non-convertible gwf cell is less than the top of each delay
    !! interbed cell.
    !!
    !<
    procedure, private :: csub_delay_head_check

    ! methods for tables
    !> @brief Initialize optional tables
    !!
    !! Subroutine to initialize optional tables. Tables include:
    !!   o delay interbeds convergence tables
    !!
    !<
    procedure, private :: csub_initialize_tables
    !
    ! -- methods for observations
    !> @brief Determine if observations are supported.
    !!
    !! Function to determine if observations are supported by the CSUB package.
    !! Observations are supported by the CSUB package.
    !!
    !<
    procedure, public :: csub_obs_supported
    !> @brief Define the observation types available in the package
    !!
    !! Method to define the observation types available in the CSUB package.
    !!
    !<
    procedure, public :: csub_df_obs
    !> @brief Read and prepare the observations
    !!
    !! Method to read and prepare the observations for the CSUB package.
    !!
    !<
    procedure, private :: csub_rp_obs
    !> @brief Set the observations for this time step
    !!
    !! Method to set the CSUB package observations for this time step.
    !!
    !<
    procedure, public :: csub_bd_obs
  end type GwfCsubType

  interface

    ! -- memory allocation and deallocation (submodules/gwf-csub-memory.f90)

    module subroutine csub_allocate_scalars(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_allocate_scalars

    module subroutine csub_allocate_arrays(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_allocate_arrays

    module subroutine csub_da(this)
      class(GwfCsubType) :: this
    end subroutine csub_da

    ! -- input processing (submodules/gwf-csub-input.f90)

    module subroutine source_options(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine source_options

    module subroutine log_options(this, warn_estress_lag)
      class(GwfCsubType), intent(inout) :: this
      logical(LGP), intent(in) :: warn_estress_lag
    end subroutine log_options

    module subroutine csub_source_dimensions(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_source_dimensions

    module subroutine csub_source_griddata(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_source_griddata

    module subroutine csub_source_packagedata(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_source_packagedata

    module subroutine csub_print_packagedata(this)
      class(GwfCsubType) :: this
    end subroutine csub_print_packagedata

    module subroutine csub_initialize_tables(this)
      class(GwfCsubType) :: this
    end subroutine csub_initialize_tables

    module subroutine define_listlabel(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine define_listlabel

    ! -- initial states and stresses (submodules/gwf-csub-stress.f90)

    module subroutine csub_cg_calc_stress(this, nodes, hnew)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: nodes !< number of active model nodes
      real(DP), dimension(nodes), intent(in) :: hnew !< current head
    end subroutine csub_cg_calc_stress

    module subroutine csub_set_initial_state(this, nodes, hnew)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: nodes !< number of active model nodes
      real(DP), dimension(nodes), intent(in) :: hnew !< current heads
    end subroutine csub_set_initial_state

    ! -- coarse-grained storage (submodules/gwf-csub-coarse.f90)

    module subroutine csub_cg_fc(this, node, tled, area, hcell, hcellold, hcof, &
                                 rhs)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: tled !< recripicol of the time step length
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head
      real(DP), intent(in) :: hcellold !< previous head
      real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
      real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    end subroutine csub_cg_fc

    module subroutine csub_cg_fn(this, node, tled, area, hcell, hcof, rhs)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: node !< node number
      real(DP), intent(in) :: tled !< reciprocal of the time step length
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
      real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    end subroutine csub_cg_fn

    module subroutine csub_cg_calc_sske(this, n, sske, hcell)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: n !< cell node number
      real(DP), intent(inout) :: sske !< coarse grained Sske
      real(DP), intent(in) :: hcell !< current head in cell
    end subroutine csub_cg_calc_sske

    module subroutine csub_cg_calc_comp(this, node, hcell, hcellold, comp)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: comp !< coarse-grained compaction
    end subroutine csub_cg_calc_comp

    module subroutine csub_cg_update(this, node)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: node !< cell node number
    end subroutine csub_cg_update

    module subroutine csub_cg_wcomp_fc(this, node, tled, area, hcell, hcellold, &
                                       hcof, rhs)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: tled !< reciprocal of the time step length
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
      real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    end subroutine csub_cg_wcomp_fc

    module subroutine csub_cg_wcomp_fn(this, node, tled, area, hcell, hcellold, &
                                       hcof, rhs)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: tled !< reciprocal of the time step length
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: hcof !< coarse-grained A matrix entry
      real(DP), intent(inout) :: rhs !< coarse-grained right-hand side entry
    end subroutine csub_cg_wcomp_fn

    ! -- no-delay interbed storage (submodules/gwf-csub-nodelay.f90)

    module subroutine csub_nodelay_update(this, i)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: i
    end subroutine csub_nodelay_update

    module subroutine csub_nodelay_fc(this, ib, hcell, hcellold, rho1, rho2, &
                                      rhs, argtled)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head in the cell
      real(DP), intent(in) :: hcellold !< previous head in the cell
      real(DP), intent(inout) :: rho1 !< current storage coefficient value using Sske
      real(DP), intent(inout) :: rho2 !< current storage coefficient value based on Ssk
      real(DP), intent(inout) :: rhs !< no-delay interbed contribution to the right-hand side
      real(DP), intent(in), optional :: argtled !< optional reciprocal of the time step length
    end subroutine csub_nodelay_fc

    module subroutine csub_nodelay_calc_comp(this, ib, hcell, hcellold, comp, &
                                             rho1, rho2)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head for the cell
      real(DP), intent(in) :: hcellold !< previous head for the cell
      real(DP), intent(inout) :: comp !< no-delay interbed compaction
      real(DP), intent(inout) :: rho1 !< current storage coefficient based on Sske
      real(DP), intent(inout) :: rho2 !< current storage coefficient based on Ssk
    end subroutine csub_nodelay_calc_comp

    module subroutine csub_interbed_fc(this, ib, node, area, hcell, hcellold, &
                                       hcof, rhs)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: hcof !< interbed A matrix entry
      real(DP), intent(inout) :: rhs !< interbed right-hand side
    end subroutine csub_interbed_fc

    module subroutine csub_interbed_fn(this, ib, node, area, hcell, hcellold, &
                                       hcof, rhs)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in a cell
      real(DP), intent(in) :: hcellold !< previous head in a cell
      real(DP), intent(inout) :: hcof !< interbed A matrix entry
      real(DP), intent(inout) :: rhs !< interbed right-hand side entry
    end subroutine csub_interbed_fn

    module subroutine csub_nodelay_wcomp_fc(this, ib, node, tled, area, hcell, &
                                            hcellold, hcof, rhs)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: tled !< reciprocal of time step length
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: hcof !< no-delay A matrix entry
      real(DP), intent(inout) :: rhs !< no-delay right-hand side entry
    end subroutine csub_nodelay_wcomp_fc

    module subroutine csub_nodelay_wcomp_fn(this, ib, node, tled, area, hcell, &
                                            hcellold, hcof, rhs)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: node !< cell node number
      real(DP), intent(in) :: tled !< reciprocal of time step length
      real(DP), intent(in) :: area !< horizontal cell area
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: hcof !< no-delay A matrix entry
      real(DP), intent(inout) :: rhs !< no-delay right-hand side entry
    end subroutine csub_nodelay_wcomp_fn

    ! -- delay interbeds (submodules/gwf-csub-delay.f90)

    module subroutine csub_delay_head_check(this, ib)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
    end subroutine csub_delay_head_check

    module subroutine csub_delay_sln(this, ib, hcell, update)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head in a cell
      logical(LGP), intent(in), optional :: update !< optional logical variable indicating
    end subroutine csub_delay_sln

    module subroutine csub_delay_init_zcell(this, ib)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
    end subroutine csub_delay_init_zcell

    module subroutine csub_delay_calc_stress(this, ib, hcell)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head in a cell
    end subroutine csub_delay_calc_stress

    module subroutine csub_delay_calc_ssksske(this, ib, n, hcell, ssk, sske, &
                                              dsskde, wfac)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: n !< delay interbed cell number
      real(DP), intent(in) :: hcell !< current head in a cell
      real(DP), intent(inout) :: ssk !< delay interbed skeletal specific storage
      real(DP), intent(inout) :: sske !< delay interbed elastic skeletal specific storage
      real(DP), intent(out), optional :: dsskde !< d(ssk)/d(effective stress)
      real(DP), intent(out), optional :: wfac !< inelastic weight (0 elastic, 1 inelastic) for the budget split
    end subroutine csub_delay_calc_ssksske

    module subroutine csub_delay_assemble(this, ib, hcell)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head in a cell
    end subroutine csub_delay_assemble

    module subroutine csub_delay_assemble_fc(this, ib, n, hcell, aii, au, al, r)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: n !< delay interbed cell number
      real(DP), intent(in) :: hcell !< current head in a cell
      real(DP), intent(inout) :: aii !< diagonal in the A matrix
      real(DP), intent(inout) :: au !< upper term in the A matrix
      real(DP), intent(inout) :: al !< lower term in the A matrix
      real(DP), intent(inout) :: r !< right-hand side term
    end subroutine csub_delay_assemble_fc

    module subroutine csub_delay_assemble_fn(this, ib, n, hcell, aii, au, al, r)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: n !< delay interbed cell number
      real(DP), intent(in) :: hcell !< current head in a cell
      real(DP), intent(inout) :: aii !< diagonal in the A matrix
      real(DP), intent(inout) :: au !< upper term in the A matrix
      real(DP), intent(inout) :: al !< lower term in the A matrix
      real(DP), intent(inout) :: r !< right-hand side term
    end subroutine csub_delay_assemble_fn

    module subroutine csub_delay_calc_sat(this, node, idelay, n, hcell, &
                                          hcellold, snnew, snold)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: node !< cell node number
      integer(I4B), intent(in) :: idelay !< delay interbed number
      integer(I4B), intent(in) :: n !< delay interbed cell number
      real(DP), intent(in) :: hcell !< current head in delay interbed cell n
      real(DP), intent(in) :: hcellold !< previous head in delay interbed cell n
      real(DP), intent(inout) :: snnew !< current saturation in delay interbed cell n
      real(DP), intent(inout) :: snold !< previous saturation in delay interbed cell n
    end subroutine csub_delay_calc_sat

    module function csub_delay_calc_sat_derivative(this, node, idelay, n, &
                                                   hcell) result(satderv)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: node !< cell node number
      integer(I4B), intent(in) :: idelay !< delay interbed number
      integer(I4B), intent(in) :: n !< delay interbed cell number
      real(DP), intent(in) :: hcell !< current head in delay interbed cell n
      real(DP) :: satderv
    end function csub_delay_calc_sat_derivative

    module subroutine csub_delay_calc_dstor(this, ib, hcell, stoe, stoi)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(inout) :: stoe !< current elastic storage change in delay interbed
      real(DP), intent(inout) :: stoi !< current inelastic storage change in delay interbed
    end subroutine csub_delay_calc_dstor

    module subroutine csub_delay_calc_wcomp(this, ib, dwc)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(inout) :: dwc !< current water compressibility change in delay interbed
    end subroutine csub_delay_calc_wcomp

    module subroutine csub_delay_calc_comp(this, ib, hcell, hcellold, comp, &
                                           compi, compe)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP), intent(in) :: hcellold !< previous head in cell
      real(DP), intent(inout) :: comp !< compaction in delay interbed
      real(DP), intent(inout) :: compi !< inelastic compaction in delay interbed
      real(DP), intent(inout) :: compe !< elastic compaction in delay interbed
    end subroutine csub_delay_calc_comp

    module subroutine csub_delay_update(this, ib)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
    end subroutine csub_delay_update

    module subroutine csub_delay_fc(this, ib, hcof, rhs)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      real(DP), intent(inout) :: hcof !< head dependent coefficient
      real(DP), intent(inout) :: rhs !< right-hand side
    end subroutine csub_delay_fc

    module function csub_calc_delay_flow(this, ib, n, hcell) result(q)
      class(GwfCsubType), intent(inout) :: this
      integer(I4B), intent(in) :: ib !< interbed number
      integer(I4B), intent(in) :: n !< delay interbed cell
      real(DP), intent(in) :: hcell !< current head in cell
      real(DP) :: q
    end function csub_calc_delay_flow

    ! -- convergence check and budget (submodules/gwf-csub-budget.f90)

    module subroutine csub_cc(this, innertot, kiter, iend, icnvgmod, nodes, &
                              hnew, hold, cpak, ipak, dpak)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: innertot !< total number of inner iterations
      integer(I4B), intent(in) :: kiter !< outer iteration number
      integer(I4B), intent(in) :: iend !< flag indicating if it is the last iteration
      integer(I4B), intent(in) :: icnvgmod !< flag indicating if the solution is considered converged
      integer(I4B), intent(in) :: nodes !< number of active nodes
      real(DP), dimension(nodes), intent(in) :: hnew !< current gwf head
      real(DP), dimension(nodes), intent(in) :: hold !< gwf for previous time step
      character(len=LENPAKLOC), intent(inout) :: cpak !< string location of the maximum change in csub package
      integer(I4B), intent(inout) :: ipak !< node with the maximum change in csub package
      real(DP), intent(inout) :: dpak !< maximum change in csub package
    end subroutine csub_cc

    module subroutine csub_cq(this, nodes, hnew, hold, isuppress_output, flowja)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: nodes !< number of active model nodes
      real(DP), intent(in), dimension(nodes) :: hnew !< current head
      real(DP), intent(in), dimension(nodes) :: hold !< head for the previous time step
      integer(I4B), intent(in) :: isuppress_output !< flag indicating if budget output should be suppressed
      real(DP), dimension(:), contiguous, intent(inout) :: flowja
    end subroutine csub_cq

    module subroutine csub_bd(this, isuppress_output, model_budget)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: isuppress_output
      type(BudgetType), intent(inout) :: model_budget !< model budget object
    end subroutine csub_bd

    module subroutine csub_save_model_flows(this, icbcfl, icbcun)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: icbcfl !< flag to output budget data
      integer(I4B), intent(in) :: icbcun !< unit number for cell-by-cell file
    end subroutine csub_save_model_flows

    ! -- output (submodules/gwf-csub-output.f90)

    module subroutine csub_fp(this)
      class(GwfCsubType) :: this
    end subroutine csub_fp

    module subroutine csub_ot_dv(this, idvfl, idvprint)
      class(GwfCsubType) :: this
      integer(I4B), intent(in) :: idvfl !< flag to save dependent variable data
      integer(I4B), intent(in) :: idvprint !< flag to print dependent variable data
    end subroutine csub_ot_dv

    module subroutine csub_cg_chk_stress(this)
      class(GwfCsubType) :: this
    end subroutine csub_cg_chk_stress

    ! -- observations (submodules/gwf-csub-obs.f90)

    module logical function csub_obs_supported(this)
      class(GwfCsubType) :: this
    end function csub_obs_supported

    module subroutine csub_df_obs(this)
      class(GwfCsubType) :: this
    end subroutine csub_df_obs

    module subroutine csub_bd_obs(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_bd_obs

    module subroutine csub_rp_obs(this)
      class(GwfCsubType), intent(inout) :: this
    end subroutine csub_rp_obs

    module subroutine csub_process_obsID(obsrv, dis, inunitobs, iout)
      type(ObserveType), intent(inout) :: obsrv !< observation type
      class(DisBaseType), intent(in) :: dis !< pointer to the model discretization
      integer(I4B), intent(in) :: inunitobs !< unit number of the observation file
      integer(I4B), intent(in) :: iout !< unit number to the model listing file
    end subroutine csub_process_obsID

  end interface

contains

  !> @brief Create a new package object
  !!
  !!  Create a new CSUB object
  !!
  !<
  subroutine csub_cr(csubobj, name_model, mempath, istounit, stoPckName, inunit, &
                     iout)
    ! -- dummy variables
    type(GwfCsubType), pointer :: csubobj !< pointer to default package type
    character(len=*), intent(in) :: name_model !< model name
    character(len=*), intent(in) :: mempath !< input context mem path
    integer(I4B), intent(in) :: inunit !< unit number of csub input file
    integer(I4B), intent(in) :: istounit !< unit number of storage package
    character(len=*), intent(in) :: stoPckName !< name of the storage package
    integer(I4B), intent(in) :: iout !< unit number of lst output file
    ! -- local variables
    !
    ! -- allocate the object and assign values to object variables
    allocate (csubobj)

    ! -- create name and memory path
    call csubobj%set_names(1, name_model, 'CSUB', 'CSUB', mempath)
    !
    ! -- Allocate scalars
    call csubobj%csub_allocate_scalars()
    !
    ! -- Create memory path to variables from STO package
    csubobj%stoMemPath = create_mem_path(name_model, stoPckName)
    !
    ! -- Set variables
    csubobj%istounit = istounit
    csubobj%inunit = inunit
    csubobj%iout = iout
  end subroutine csub_cr

  !> @brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the CSUB package.
  !!
  !<
  subroutine csub_ar(this, dis, ibound)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use ConstantsModule, only: LINELENGTH
    use KindModule, only: I4B
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    class(DisBaseType), pointer, intent(in) :: dis !< model discretization
    integer(I4B), dimension(:), pointer, contiguous :: ibound !< model ibound array
    ! -- local variables
    character(len=20) :: cellid
    integer(I4B) :: idelay
    integer(I4B) :: ib
    integer(I4B) :: node
    integer(I4B) :: n
    integer(I4B) :: istoerr
    real(DP) :: top
    real(DP) :: bot
    real(DP) :: thick
    real(DP) :: cg_ske_cr
    real(DP) :: theta
    real(DP) :: v
    real(DP) :: vtot
    real(DP) :: cell_thickness
    real(DP) :: overshoot
    real(DP) :: f
    real(DP) :: rval
    ! -- format
    character(len=*), parameter :: fmtcsub = &
      "(1x,/1x,'CSUB -- COMPACTION PACKAGE, VERSION 1.15, 7/27/2026', &
     &' INPUT READ FROM MEMPATH: ', A, /)"
    !
    ! --print a message identifying the csub package.
    write (this%iout, fmtcsub) this%input_mempath
    !
    ! -- store pointers to arguments that were passed in
    this%dis => dis
    this%ibound => ibound
    !
    ! -- create obs package
    call obs_cr(this%obs, this%inobspkg)
    !
    ! -- source csub options
    call this%source_options()
    !
    ! -- source the csub dimensions
    call this%source_dimensions()
    !
    ! - observation data
    call this%obs%obs_ar()
    !
    ! -- terminate if errors dimensions block data
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    ! -- Allocate arrays in
    call this%csub_allocate_arrays()
    !
    ! -- source griddata
    call this%csub_source_griddata()
    !
    ! -- evaluate the coarse-grained material properties and if
    !    non-zero specific storage values are specified in the
    !    STO package
    istoerr = 0
    do node = 1, this%dis%nodes
      call this%dis%noder_to_string(node, cellid)
      cg_ske_cr = this%cg_ske_cr(node)
      theta = this%cg_thetaini(node)
      !
      ! -- coarse-grained storage error condition
      if (cg_ske_cr < DZERO) then
        write (errmsg, '(a,g0,a,1x,a,1x,a,a)') &
          'Coarse-grained material CG_SKE_CR (', cg_ske_cr, ') is less', &
          'than zero in cell', trim(adjustl(cellid)), '.'
      end if
      !
      ! -- storage (STO) package error condition
      if (this%stoss(node) /= DZERO) then
        istoerr = 1
      end if
      !
      ! -- porosity error condition
      if (theta > DONE .or. theta < DZERO) then
        write (errmsg, '(a,g0,a,1x,a,1x,a,a)') &
          'Coarse-grained material THETA (', theta, ') is less', &
          'than zero or greater than 1 in cell', trim(adjustl(cellid)), '.'
      end if
    end do
    !
    ! -- write single message if storage (STO) package has non-zero specific
    !    storage values
    if (istoerr /= 0) then
      write (errmsg, '(a,3(1x,a))') &
        'Specific storage values in the storage (STO) package must', &
        'be zero in all active cells when using the', &
        trim(adjustl(this%packName)), &
        'package.'
      call store_error(errmsg)
    end if
    !
    ! -- source interbed data
    if (this%ninterbeds > 0) then
      call this%csub_source_packagedata()
    end if
    !
    ! setup package convergence tables
    call this%csub_initialize_tables()
    !
    ! -- calculate the coarse-grained material thickness without the interbeds
    do node = 1, this%dis%nodes
      top = this%dis%top(node)
      bot = this%dis%bot(node)
      this%cg_thickini(node) = top - bot
      this%cell_thick(node) = top - bot
    end do
    !
    ! -- subtract the interbed thickness from aquifer thickness
    do ib = 1, this%ninterbeds
      node = this%nodelist(ib)
      idelay = this%idelay(ib)
      if (idelay == 0) then
        v = this%thickini(ib)
      else
        v = this%rnb(ib) * this%thickini(ib)
      end if
      this%cg_thickini(node) = this%cg_thickini(node) - v
    end do
    !
    ! -- evaluate if the interbed thicknesses in a cell exceed the cell
    !    thickness (cg_thickini < 0). A small excess (numerical roundoff in the
    !    summed interbed thicknesses) is resolved by proportionally scaling the
    !    cell's interbed thicknesses to fit; a larger excess is a genuine
    !    over-specification and is an error
    do node = 1, this%dis%nodes
      thick = this%cg_thickini(node)
      if (thick >= DZERO) cycle
      cell_thickness = this%cell_thick(node)
      overshoot = -thick
      if (overshoot <= DEM5 * cell_thickness) then
        ! -- roundoff: scale the cell's interbed thicknesses so they fill the
        !    cell exactly (cell_thickness - thick is the summed interbed
        !    thickness, which exceeds cell_thickness by overshoot)
        f = cell_thickness / (cell_thickness - thick)
        do ib = 1, this%ninterbeds
          if (node /= this%nodelist(ib)) cycle
          this%thickini(ib) = this%thickini(ib) * f
          if (this%iupdatematprop /= 0) then
            this%thick(ib) = this%thick(ib) * f
          end if
          idelay = this%idelay(ib)
          if (idelay /= 0) then
            rval = this%thickini(ib) / real(this%ndelaycells, DP)
            do n = 1, this%ndelaycells
              this%dbdzini(n, idelay) = rval
              this%dbdz(n, idelay) = rval
              this%dbdz0(n, idelay) = rval
            end do
            ! -- recompute delay-bed cell elevations from the scaled thickness
            call this%csub_delay_init_zcell(ib)
          end if
        end do
        this%cg_thickini(node) = DZERO
      else
        call this%dis%noder_to_string(node, cellid)
        write (errmsg, '(a,g0,a,1x,a,a)') &
          'Coarse grained material thickness is less than zero (', &
          thick, ') in cell', trim(adjustl(cellid)), '. Interbed thicknesses:'

        vtot = DZERO
        do ib = 1, this%ninterbeds
          if (node /= this%nodelist(ib)) then
            cycle
          end if
          idelay = this%idelay(ib)
          v = this%thickini(ib)
          if (idelay /= 0) then
            v = v * this%rnb(ib)
          end if
          vtot = vtot + v
          write (errmsg, '(a,1x,a,i0,a,g0)') &
            trim(adjustl(errmsg)), &
            'icbno(', ib, ')=', v
        end do
        write (errmsg, '(a,a,g0,a)') &
          trim(adjustl(errmsg)), &
          '. Total interbed thickness=', vtot, '.'
        call store_error(errmsg)
      end if
    end do
    !
    ! -- terminate if errors griddata, packagedata blocks, TDIS, or STO data
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- set current coarse-grained thickness (cg_thick) and
    !    current coarse-grained porosity (cg_theta). Only needed
    !    if updating material properties
    if (this%iupdatematprop /= 0) then
      do node = 1, this%dis%nodes
        this%cg_thick(node) = this%cg_thickini(node)
        this%cg_theta(node) = this%cg_thetaini(node)
      end do
    end if
  end subroutine csub_ar

  !> @brief Read and prepare stress period data for package
  !!
  !!  Method reads and prepares stress period data for the CSUB package.
  !!  The overlying geostatic stress (sig0) is the only stress period data
  !!  read by the CSUB package.
  !!
  !<
  subroutine csub_rp(this)
    ! -- modules
    use TdisModule, only: kper
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    ! -- local variables
    integer(I4B), dimension(:, :), pointer, contiguous :: cellids
    integer(I4B), dimension(:), pointer, contiguous :: cellid
    integer(I4B), pointer :: iper
    integer(I4B) :: n, nodeu, noder
    character(len=LINELENGTH) :: title, text
    character(len=20) :: cellstr
    logical(LGP) :: found
    ! -- formats
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',a,'S FROM LAST STRESS PERIOD')"

    call mem_setptr(iper, 'IPER', this%input_mempath)
    if (iper /= kper) then
      write (this%iout, fmtlsp) trim(this%filtyp)
      call this%csub_rp_obs()
      return
    end if

    call mem_setptr(cellids, 'CELLID', this%input_mempath)
    call mem_set_value(this%nbound, 'NBOUND', this%input_mempath, &
                       found, release=.false.)

    ! -- setup table for period data
    if (this%iprpak /= 0) then
      ! -- reset the input table object
      title = 'CSUB'//' PACKAGE ('// &
              trim(adjustl(this%packName))//') DATA FOR PERIOD'
      write (title, '(a,1x,i6)') trim(adjustl(title)), kper
      call table_cr(this%inputtab, this%packName, title)
      call this%inputtab%table_df(1, 2, this%iout, finalize=.FALSE.)
      text = 'CELLID'
      call this%inputtab%initialize_column(text, 20)
      text = 'SIG0'
      call this%inputtab%initialize_column(text, 15, alignment=TABLEFT)
    end if

    ! -- update nodelist
    do n = 1, this%nbound

      ! -- set cellid
      cellid => cellids(:, n)

      ! -- set user node number
      if (this%dis%ndim == 1) then
        nodeu = cellid(1)
      elseif (this%dis%ndim == 2) then
        nodeu = get_node(cellid(1), 1, cellid(2), &
                         this%dis%mshape(1), 1, &
                         this%dis%mshape(2))
      else
        nodeu = get_node(cellid(1), cellid(2), cellid(3), &
                         this%dis%mshape(1), &
                         this%dis%mshape(2), &
                         this%dis%mshape(3))
      end if

      ! -- set noder
      noder = this%dis%get_nodenumber(nodeu, 1)
      if (noder <= 0) then
        cycle
      end if

      this%nodelistsig0(n) = noder

      ! -- write line to table
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(noder, cellstr)
        call this%inputtab%add_term(cellstr)
        call this%inputtab%add_term(this%sig0(n))
      end if
    end do
    !
    ! -- terminate if errors encountered
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
    !
    ! -- finalize the table
    if (this%iprpak /= 0) then
      call this%inputtab%finalize_table()
    end if
    !
    ! -- read observations
    call this%csub_rp_obs()
  end subroutine csub_rp

  !> @brief Advance the package
  !!
  !!  Advance data in the CSUB package. The method sets data for the previous
  !!  time step to the current value for the data (e.g., HOLD = HNEW). The
  !!  method also calls the method to initialize the initial stress conditions
  !!  if this is the first transient stress period.
  !!
  !<
  subroutine csub_ad(this, nodes, hnew)
    ! -- modules
    use TdisModule, only: nper, kper
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: nodes !< number of active model nodes
    real(DP), dimension(nodes), intent(in) :: hnew !< current head
    ! -- local variables
    integer(I4B) :: ib
    integer(I4B) :: n
    integer(I4B) :: idelay
    integer(I4B) :: node
    real(DP) :: h
    real(DP) :: es
    real(DP) :: pcs
    !
    ! -- evaluate if steady-state stress periods are specified for more
    !    than the first and last stress period if interbeds are simulated
    if (this%ninterbeds > 0) then
      if (kper > 1 .and. kper < nper) then
        if (this%gwfiss /= 0) then
          write (errmsg, '(a,i0,a,1x,a,1x,a,1x,i0,1x,a)') &
            'Only the first and last (', nper, ')', &
            'stress period can be steady if interbeds are simulated.', &
            'Stress period', kper, 'has been defined to be steady state.'
          call store_error(errmsg, terminate=.TRUE.)
        end if
      end if
    end if
    !
    ! -- set initial states
    if (this%initialized == 0) then
      if (this%gwfiss == 0) then
        call this%csub_set_initial_state(nodes, hnew)
      end if
    end if
    !
    ! -- update state variables
    !
    ! -- coarse-grained materials
    do node = 1, nodes
      this%cg_comp(node) = DZERO
      this%cg_es0(node) = this%cg_es(node)
      if (this%iupdatematprop /= 0) then
        this%cg_thick0(node) = this%cg_thick(node)
        this%cg_theta0(node) = this%cg_theta(node)
      end if
    end do
    !
    ! -- interbeds
    do ib = 1, this%ninterbeds
      idelay = this%idelay(ib)
      !
      ! -- update common terms for no-delay and delay interbeds
      this%comp(ib) = DZERO
      node = this%nodelist(ib)
      if (this%initialized /= 0) then
        es = this%cg_es(node)
        pcs = this%pcs(ib)
        if (es > pcs) then
          this%pcs(ib) = es
        end if
      end if
      if (this%iupdatematprop /= 0) then
        this%thick0(ib) = this%thick(ib)
        this%theta0(ib) = this%theta(ib)
      end if
      !
      ! -- update delay interbed terms
      if (idelay /= 0) then
        !
        ! -- update state if previous period was steady state
        if (kper > 1) then
          if (this%gwfiss0 /= 0) then
            node = this%nodelist(ib)
            h = hnew(node)
            do n = 1, this%ndelaycells
              this%dbh(n, idelay) = h
            end do
          end if
        end if
        !
        ! -- update preconsolidation stress, stresses, head, dbdz0, and theta0
        do n = 1, this%ndelaycells
          ! update preconsolidation stress
          if (this%initialized /= 0) then
            if (this%dbes(n, idelay) > this%dbpcs(n, idelay)) then
              this%dbpcs(n, idelay) = this%dbes(n, idelay)
            end if
          end if
          this%dbh0(n, idelay) = this%dbh(n, idelay)
          this%dbes0(n, idelay) = this%dbes(n, idelay)
          if (this%iupdatematprop /= 0) then
            this%dbdz0(n, idelay) = this%dbdz(n, idelay)
            this%dbtheta0(n, idelay) = this%dbtheta(n, idelay)
          end if
        end do
      end if
    end do
    !
    ! -- set gwfiss0
    this%gwfiss0 = this%gwfiss
    !
    ! -- For each observation, push simulated value and corresponding
    !    simulation time from "current" to "preceding" and reset
    !    "current" value.
    call this%obs%obs_ad()
  end subroutine csub_ad

  !> @brief Fill A and r for the package
  !!
  !!  Fill the coefficient matrix and right-hand side with the CSUB package terms.
  !!
  !<
  subroutine csub_fc(this, kiter, hold, hnew, matrix_sln, idxglo, rhs)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: kiter !< outer iteration numbed
    real(DP), intent(in), dimension(:) :: hold !< previous heads
    real(DP), intent(in), dimension(:) :: hnew !< current heads
    class(MatrixBaseType), pointer :: matrix_sln !< A matrix
    integer(I4B), intent(in), dimension(:) :: idxglo !< global index model to solution
    real(DP), intent(inout), dimension(:) :: rhs !< right-hand side
    ! -- local variables
    integer(I4B) :: ib
    integer(I4B) :: node
    integer(I4B) :: idiag
    integer(I4B) :: idelay
    real(DP) :: tled
    real(DP) :: area
    real(DP) :: hcof
    real(DP) :: rhsterm
    real(DP) :: comp
    !
    ! -- update geostatic load calculation
    call this%csub_cg_calc_stress(this%dis%nodes, hnew)
    !
    ! -- formulate csub terms
    if (this%gwfiss == 0) then
      !
      ! -- initialize tled
      tled = DONE / delt
      !
      ! -- coarse-grained storage
      do node = 1, this%dis%nodes
        idiag = this%dis%con%ia(node)
        area = this%dis%get_area(node)
        !
        ! -- skip inactive cells
        if (this%ibound(node) < 1) cycle
        !
        ! -- update coarse-grained material properties
        if (this%iupdatematprop /= 0) then
          if (this%ieslag == 0) then
            !
            ! -- calculate compaction
            call this%csub_cg_calc_comp(node, hnew(node), hold(node), comp)
            this%cg_comp(node) = comp
            !
            ! -- update coarse-grained thickness and void ratio
            call this%csub_cg_update(node)
          end if
        end if
        !
        ! -- calculate coarse-grained storage terms
        call this%csub_cg_fc(node, tled, area, hnew(node), hold(node), &
                             hcof, rhsterm)
        !
        ! -- add coarse-grained storage terms to amat and rhs for coarse-grained storage
        call matrix_sln%add_value_pos(idxglo(idiag), hcof)
        rhs(node) = rhs(node) + rhsterm
        !
        ! -- calculate coarse-grained water compressibility
        !    storage terms
        if (this%brg /= DZERO) then
          call this%csub_cg_wcomp_fc(node, tled, area, hnew(node), hold(node), &
                                     hcof, rhsterm)
          !
          ! -- add water compression storage terms to amat and rhs for
          !   coarse-grained storage
          call matrix_sln%add_value_pos(idxglo(idiag), hcof)
          rhs(node) = rhs(node) + rhsterm
        end if
      end do
      !
      ! -- interbed storage
      if (this%ninterbeds /= 0) then
        !
        ! -- calculate the contribution of interbeds to the
        !    groundwater flow equation
        do ib = 1, this%ninterbeds
          node = this%nodelist(ib)
          idelay = this%idelay(ib)
          idiag = this%dis%con%ia(node)
          area = this%dis%get_area(node)
          call this%csub_interbed_fc(ib, node, area, hnew(node), hold(node), &
                                     hcof, rhsterm)
          call matrix_sln%add_value_pos(idxglo(idiag), hcof)
          rhs(node) = rhs(node) + rhsterm
          !
          ! -- calculate interbed water compressibility terms
          if (.not. is_close(this%brg, DZERO) .and. idelay == 0) then
            call this%csub_nodelay_wcomp_fc(ib, node, tled, area, &
                                            hnew(node), hold(node), &
                                            hcof, rhsterm)
            !
            ! -- add water compression storage terms to amat and rhs for interbed
            call matrix_sln%add_value_pos(idxglo(idiag), hcof)
            rhs(node) = rhs(node) + rhsterm
          end if
        end do
      end if
    end if
    !
    ! -- terminate if errors encountered when updating material properties
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine csub_fc

  !> @brief Fill Newton-Raphson terms in A and r for the package
  !!
  !!  Fill the coefficient matrix and right-hand side with CSUB package
  !!  with Newton-Raphson terms.
  !<
  subroutine csub_fn(this, kiter, hold, hnew, matrix_sln, idxglo, rhs)
    ! -- modules
    use TdisModule, only: delt
    ! -- dummy variables
    class(GwfCsubType) :: this
    integer(I4B), intent(in) :: kiter !< outer iteration number
    real(DP), intent(in), dimension(:) :: hold !< previous heads
    real(DP), intent(in), dimension(:) :: hnew !< current heads
    class(MatrixBaseType), pointer :: matrix_sln !< A matrix
    integer(I4B), intent(in), dimension(:) :: idxglo !< global index model to solution
    real(DP), intent(inout), dimension(:) :: rhs !< right-hand side
    ! -- local variables
    integer(I4B) :: idelay
    integer(I4B) :: node
    integer(I4B) :: idiag
    integer(I4B) :: ib
    real(DP) :: tled
    real(DP) :: area
    real(DP) :: hcof
    real(DP) :: rhsterm
    !
    ! -- formulate csub terms
    if (this%gwfiss == 0) then
      tled = DONE / delt
      !
      ! -- coarse-grained storage
      do node = 1, this%dis%nodes
        idiag = this%dis%con%ia(node)
        area = this%dis%get_area(node)
        !
        ! -- skip inactive cells
        if (this%ibound(node) < 1) cycle
        !
        ! -- calculate coarse-grained storage newton terms
        call this%csub_cg_fn(node, tled, area, &
                             hnew(node), hcof, rhsterm)
        !
        ! -- add coarse-grained storage newton terms to amat and rhs for
        !   coarse-grained storage
        call matrix_sln%add_value_pos(idxglo(idiag), hcof)
        rhs(node) = rhs(node) + rhsterm
        !
        ! -- calculate coarse-grained water compressibility storage
        !    newton terms
        if (this%brg /= DZERO) then
          call this%csub_cg_wcomp_fn(node, tled, area, hnew(node), hold(node), &
                                     hcof, rhsterm)
          !
          ! -- add water compression storage newton terms to amat and rhs for
          !    coarse-grained storage
          call matrix_sln%add_value_pos(idxglo(idiag), hcof)
          rhs(node) = rhs(node) + rhsterm
        end if
      end do
      !
      ! -- interbed storage
      if (this%ninterbeds /= 0) then
        !
        ! -- calculate the interbed newton terms for the
        !    groundwater flow equation
        do ib = 1, this%ninterbeds
          idelay = this%idelay(ib)
          node = this%nodelist(ib)
          !
          ! -- skip inactive cells
          if (this%ibound(node) < 1) cycle
          !
          ! -- calculate interbed newton terms
          idiag = this%dis%con%ia(node)
          area = this%dis%get_area(node)
          call this%csub_interbed_fn(ib, node, area, hnew(node), hold(node), &
                                     hcof, rhsterm)
          !
          ! -- add interbed newton terms to amat and rhs
          call matrix_sln%add_value_pos(idxglo(idiag), hcof)
          rhs(node) = rhs(node) + rhsterm
          !
          ! -- calculate interbed water compressibility terms
          if (this%brg /= DZERO .and. idelay == 0) then
            call this%csub_nodelay_wcomp_fn(ib, node, tled, area, &
                                            hnew(node), hold(node), &
                                            hcof, rhsterm)
            !
            ! -- add interbed water compression newton terms to amat and rhs
            call matrix_sln%add_value_pos(idxglo(idiag), hcof)
            rhs(node) = rhs(node) + rhsterm
          end if
        end do
      end if
    end if
  end subroutine csub_fn

  !> @brief Calculate the interbed thickness
  !!
  !! Function to calculate the interbed thickness.
  !!
  !! @return      thick               interbed thickness
  !<
  function csub_calc_interbed_thickness(this, ib) result(thick)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: ib !< interbed number
    ! -- local variables
    integer(I4B) :: idelay
    real(DP) :: thick
    !
    ! -- calculate interbed thickness
    idelay = this%idelay(ib)
    thick = this%thick(ib)
    if (idelay /= 0) then
      thick = thick * this%rnb(ib)
    end if
  end function csub_calc_interbed_thickness

  !> @brief Calculate the cell node
  !!
  !! Function to calculate elevation of the node between the specified corrected
  !! elevation zbar and the bottom elevation. If zbar is greater than the top
  !! elevation, the node elevation is halfway between the top and bottom
  !! elevations. The corrected elevation (zbar) is always greater than or
  !! equal to bottom.
  !!
  !! @return      znode               node elevation
  !<
  function csub_calc_znode(this, top, bottom, zbar) result(znode)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    real(DP), intent(in) :: top !< top of cell
    real(DP), intent(in) :: bottom !< bottom of cell
    real(DP), intent(in) :: zbar !< corrected elevation
    ! -- local variables
    real(DP) :: znode
    real(DP) :: v
    !
    ! -- calculate the node elevation
    if (zbar > top) then
      v = top
    else
      v = zbar
    end if
    znode = DHALF * (v + bottom)
  end function csub_calc_znode

  !> @brief Calculate cell saturation
  !!
  !! Method to calculate the cell saturation for the current and
  !! previous time step.
  !<
  subroutine csub_calc_sat(this, node, hcell, hcellold, snnew, snold)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: hcell !< current head
    real(DP), intent(in) :: hcellold !< previous head
    real(DP), intent(inout) :: snnew !< current saturation
    real(DP), intent(inout) :: snold !< previous saturation
    ! -- local variables
    real(DP) :: top
    real(DP) :: bot
    !
    ! -- calculate cell saturation
    if (this%stoiconv(node) /= 0) then
      top = this%dis%top(node)
      bot = this%dis%bot(node)
      snnew = sQuadraticSaturation(top, bot, hcell, this%satomega)
      snold = sQuadraticSaturation(top, bot, hcellold, this%satomega)
    else
      snnew = DONE
      snold = DONE
    end if
    if (this%ieslag /= 0) then
      snold = snnew
    end if
  end subroutine csub_calc_sat

  !> @brief Calculate the saturation derivative
  !!
  !! Function to calculate the derivative of the saturation with
  !! respect to the current head.
  !!
  !! @return      satderv              derivative of saturation
  !<
  function csub_calc_sat_derivative(this, node, hcell) result(satderv)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: hcell !< current head
    ! -- local variables
    real(DP) :: satderv
    real(DP) :: top
    real(DP) :: bot

    if (this%stoiconv(node) /= 0) then
      top = this%dis%top(node)
      bot = this%dis%bot(node)
      satderv = sQuadraticSaturationDerivative(top, bot, hcell, this%satomega)
    else
      satderv = DZERO
    end if
  end function csub_calc_sat_derivative

  !> @brief Calculate specific storage coefficient factor
  !!
  !! Method to calculate the factor that is used to calculate skeletal
  !! specific storage coefficients. Can be used for coarse-grained
  !! materials and interbeds.
  !<
  subroutine csub_calc_sfacts(this, node, bot, znode, theta, es, es0, geo, fact)
    ! -- dummy variables
    class(GwfCsubType), intent(inout) :: this
    integer(I4B), intent(in) :: node !< cell node number
    real(DP), intent(in) :: bot !
    real(DP), intent(in) :: znode
    real(DP), intent(in) :: theta !< porosity
    real(DP), intent(in) :: es !< current effective stress
    real(DP), intent(in) :: es0 !< previous effective stress
    real(DP), intent(in) :: geo !< geostatic stress (regularization reference)
    real(DP), intent(inout) :: fact !< skeletal storage coefficient factor (1/((1+void_ratio)*bar(es)))
    ! -- local variables
    real(DP) :: esv
    real(DP) :: adjes
    real(DP) :: esfloor
    !
    ! -- initialize variables
    fact = DZERO
    if (this%ieslag /= 0) then
      esv = es0
    else
      esv = es
    end if
    !
    ! -- effective stress adjusted to the vertical node position
    adjes = adjust_eff_stress(esv, bot, znode, this%sgs(node))
    !
    ! -- smoothly floor the adjusted effective stress at stressfloor * geo so the
    !    storage factor (1/es) stays bounded and positive as es approaches zero;
    !    unchanged above the floor, disabled by STRICT_EFFECTIVE_STRESS
    if (this%istrict_stress == 0 .and. geo > DZERO) then
      esfloor = stressfloor * geo
      adjes = sQuadratic0sp(adjes, esfloor, esfloor)
    end if
    !
    ! -- calculate storage factors for the effective stress case
    fact = storage_factor(adjes, theta)
  end subroutine csub_calc_sfacts

end module GwfCsubModule
