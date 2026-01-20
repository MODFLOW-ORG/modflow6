!> @brief This module contains simulation constants
!!
!! This module contains simulation constants that are available to all
!! other modules. The variables in this module are defined at run time.
!! The module does not have any dependencies on models, exchanges, or
!! solutions in a simulation.
!!
!<
module ConstantsModule
  use KindModule
  public
  ! -- constants
  integer(I4B), parameter :: IUSERFORMATSTRIP = -99 !< default user format strip
  integer(I4B), parameter :: IUSERFORMATWRAP = 99 !< default user format wrap
  integer(I4B), parameter :: LENBIGLINE = 5000 !< maximum length of a big line
  integer(I4B), parameter :: LENHUGELINE = 50000 !< maximum length of a huge line
  integer(I4B), parameter :: LENVARNAME = 16 !< maximum length of a variable name
  integer(I4B), parameter :: LENCOMPONENTNAME = 16 !< maximum length of a component name
  integer(I4B), parameter :: LENCONTEXTNAME = 16 !< maximum length of a memory manager context
  integer(I4B), parameter :: LENANAME = 24 !< maximum length of the header text for an array
  integer(I4B), parameter :: LENSOLUTIONNAME = LENCOMPONENTNAME !< maximum length of the solution name
  integer(I4B), parameter :: LENMODELNAME = LENCOMPONENTNAME !< maximum length of the model name
  integer(I4B), parameter :: LENPACKAGENAME = LENCOMPONENTNAME !< maximum length of the package name
  integer(I4B), parameter :: LENEXCHANGENAME = LENCOMPONENTNAME !< maximum length of the exchange name
  integer(I4B), parameter :: LENBUDROWLABEL = 2 * LENPACKAGENAME + 1 !< maximum length of the rowlabel string used in the budget table
  integer(I4B), parameter :: LENMEMSEPARATOR = 1 !< maximum length of the memory path separator used, currently a '/'
  integer(I4B), parameter :: LENMEMPATH = &
                             LENCONTEXTNAME + &
                             2 * LENCOMPONENTNAME + &
                             2 * LENMEMSEPARATOR !< maximum length of the memory path
  integer(I4B), parameter :: LENMEMADDRESS = &
                             LENMEMPATH + &
                             LENMEMSEPARATOR + &
                             LENVARNAME !< maximum length of the full memory address, including variable name
  integer(I4B), parameter :: LENAUXNAME = 16 !< maximum length of a aux variable
  integer(I4B), parameter :: LENBOUNDNAME = 40 !< maximum length of a bound name
  integer(I4B), parameter :: LENBUDTXT = 16 !< maximum length of a budget component names
  integer(I4B), parameter :: LENPACKAGETYPE = 7 !< maximum length of a package type (DIS6, SFR6, CSUB6, etc.)
  integer(I4B), parameter :: LENFTYPE = 5 !< maximum length of a package type (DIS, WEL, OC, etc.)
  integer(I4B), parameter :: LENOBSNAME = 40 !< maximum length of a observation name
  integer(I4B), parameter :: LENOBSTYPE = 30 !< maximum length of a observation type (CONTINUOUS)
  integer(I4B), parameter :: LENTIMESERIESNAME = LENOBSNAME !< maximum length of a time series name
  integer(I4B), parameter :: LENTIMESERIESTEXT = 16 !< maximum length of a time series text
  integer(I4B), parameter :: LENDATETIME = 30 !< maximum length of a date time string
  integer(I4B), parameter :: LINELENGTH = 300 !< maximum length of a standard line
  integer(I4B), parameter :: LENLISTLABEL = 500 !< maximum length of a llist label
  integer(I4B), parameter :: MAXCHARLEN = max(1000, LENBIGLINE) !< maximum length of char string
  integer(I4B), parameter :: MAXOBSTYPES = 100 !< maximum number of observation types
  integer(I4B), parameter :: NAMEDBOUNDFLAG = -9 !< named bound flag
  integer(I4B), parameter :: LENPAKLOC = 34 !< maximum length of a package location
  integer(I4B), parameter :: IZERO = 0 !< integer constant zero
  integer(I4B), parameter :: IWETLAKE = 10000 !< integer constant for a dry lake
  integer(I4B), parameter :: MAXADPIT = 100 !< maximum advanced package Newton-Raphson iterations
  !
  ! -- file constants
  integer(I4B), parameter :: IUOC = 999 !< open/close file unit number
  integer(I4B), parameter :: IUSTART = 1000 !< starting file unit number
  integer(I4B), parameter :: IULAST = 10000 !< maximum file unit number (this allows for 9000 open files)
  !
  ! -- memory manager constants
  integer(I4B), public, parameter :: MAXMEMRANK = 3 !< maximum memory manager length (up to 3-dimensional arrays)
  integer(I4B), public, parameter :: LENMEMTYPE = 50 !< maximum length of a memory manager type
  !
  ! -- real constants
  real(DP), parameter :: DZERO = 0.0_DP !< real constant zero
  real(DP), parameter :: DQUARTER = 0.25_DP !< real constant 1/3
  real(DP), parameter :: DONETHIRD = 1.0_DP / 3.0_DP !< real constant 1/3
  real(DP), parameter :: DHALF = 0.5_DP !< real constant 1/2
  real(DP), parameter :: DP6 = 0.6_DP !< real constant 3/5
  real(DP), parameter :: DTWOTHIRDS = 2.0_DP / 3.0_DP !< real constant 2/3
  real(DP), parameter :: DP7 = 0.7_DP !< real constant 7/10
  real(DP), parameter :: DP9 = 0.9_DP !< real constant 9/10
  real(DP), parameter :: DP99 = 0.99_DP !< real constant 99/100
  real(DP), parameter :: DP999 = 0.999_DP !< real constant 999/1000

  real(DP), parameter :: DONE = 1.0_DP !< real constant 1
  real(DP), parameter :: D1P1 = 1.1_DP !< real constant 1.1
  real(DP), parameter :: DFIVETHIRDS = 5.0_DP / 3.0_DP !< real constant 5/3
  real(DP), parameter :: DTWO = 2.0_DP !< real constant 2
  real(DP), parameter :: DTHREE = 3.0_DP !< real constant 3
  real(DP), parameter :: DFOUR = 4.0_DP !< real constant 4
  real(DP), parameter :: DSIX = 6.0_DP !< real constant 6
  real(DP), parameter :: DEIGHT = 8.0_DP !< real constant 8
  real(DP), parameter :: DTEN = 1.0e1_DP !< real constant 10
  real(DP), parameter :: DSIXTY = 6.0e1_DP !< real constant 60
  real(DP), parameter :: DHUNDRED = 1.0e2_DP !< real constant 100

  real(DP), parameter :: DEP3 = 1.0e3_DP !< real constant 1000
  real(DP), parameter :: DEP6 = 1.0e6_DP !< real constant 1000000
  real(DP), parameter :: DEP9 = 1.0e9_DP !< real constant 1e9
  real(DP), parameter :: DEP20 = 1.0e20_DP !< real constant 1e20

  real(DP), parameter :: DHNOFLO = 1.e30_DP !< real no flow constant
  real(DP), parameter :: DHDRY = -1.e30_DP !< real dry cell constant
  real(DP), parameter :: DNODATA = 3.0e30_DP !< real no data constant

  real(DP), parameter :: DSECPERHR = 3.6e3_DP !< real constant representing number of seconds per hour (used in tdis)
  real(DP), parameter :: DHRPERDAY = 2.4e1_DP !< real constant representing number of hours per day (used in tdis)
  real(DP), parameter :: DDYPERYR = 3.6525e2_DP !< real constant representing the average number of days per year (used in tdis)
  real(DP), parameter :: DSECPERDY = 8.64e4_DP !< real constant representing the number of seconds per day (used in tdis)
  real(DP), parameter :: DSECPERYR = 3.1557600e7_DP !< real constant representing the average number of seconds per year (used in tdis)

  real(DP), parameter :: DEM1 = 1.0e-1_DP !< real constant 1e-1
  real(DP), parameter :: D5EM2 = 5.0e-2_DP !< real constant 5e-2
  real(DP), parameter :: DEM2 = 1.0e-2_DP !< real constant 1e-2
  real(DP), parameter :: DEM3 = 1.0e-3_DP !< real constant 1e-3
  real(DP), parameter :: DEM4 = 1.0e-4_DP !< real constant 1e-4
  real(DP), parameter :: DEM5 = 1.0e-5_DP !< real constant 1e-5
  real(DP), parameter :: DEM6 = 1.0e-6_DP !< real constant 1e-6
  real(DP), parameter :: DEM7 = 1.0e-7_DP !< real constant 1e-7
  real(DP), parameter :: DEM8 = 1.0e-8_DP !< real constant 1e-8
  real(DP), parameter :: DEM9 = 1.0e-9_DP !< real constant 1e-9
  real(DP), parameter :: DEM10 = 1.0e-10_DP !< real constant 1e-10
  real(DP), parameter :: DEM12 = 1.0e-12_DP !< real constant 1e-12
  real(DP), parameter :: DEM14 = 1.0e-14_DP !< real constant 1e-14
  real(DP), parameter :: DEM15 = 1.0e-15_DP !< real constant 1e-15
  real(DP), parameter :: DEM20 = 1.0e-20_DP !< real constant 1e-20
  real(DP), parameter :: DEM30 = 1.0e-30_DP !< real constant 1e-30

  real(DP), parameter :: DPREC = EPSILON(1.0_DP) !< real constant machine precision
  real(DP), parameter :: DPRECSQRT = SQRT(DPREC)
  real(DP), parameter :: DSAME = DHUNDRED * DPREC !< real constant for values that are considered
                                                  !! the same based on machine precision

  real(DP), parameter :: DLNLOW = 0.995_DP !< real constant low ratio used to calculate log mean of K
  real(DP), parameter :: DLNHIGH = 1.005_DP !< real constant high ratio used to calculate log mean of K

  real(DP), parameter :: DPI = DFOUR * ATAN(DONE) !< real constant \f$\pi\f$
  real(DP), parameter :: DTWOPI = DTWO * DFOUR * ATAN(DONE) !< real constant \f$2 \pi\f$
  real(DP), parameter :: DPIO180 = datan(DONE) / 4.5d1 !< real constant \f$\pi/180\f$

  real(DP), parameter :: DGRAVITY = 9.80665_DP !< real constant gravitational acceleration (m/(s s))
  real(DP), parameter :: DCD = 0.61_DP !< real constant weir coefficient in SI units

  character(len=10), dimension(3, 3), parameter :: &
    cidxnames = reshape( &
    ['      NODE', '          ', '          ', &
     '     LAYER', '    CELL2D', '          ', &
     '     LAYER', '       ROW', '       COL'], [3, 3]) !< cellid labels for DIS, DISV, and DISU discretizations

  !> @brief enumerator used with TimeSeriesType
  !<
  ENUM, BIND(C)
    ENUMERATOR :: UNDEFINED = 0 !< undefined interpolation
    ENUMERATOR :: STEPWISE = 1 !< stepwise interpolation
    ENUMERATOR :: LINEAR = 2 !< linear interpolation
    ENUMERATOR :: LINEAREND = 3 !< linear end interpolation
  END ENUM

  !> @brief enumerator associated with Discretization types
  !<
  ENUM, BIND(C)
    ENUMERATOR :: DISUNDEF = 0 !< undefined discretization
    ! -- 3D (base) discretizations
    ENUMERATOR :: DIS = 1 !< DIS6 discretization
    ENUMERATOR :: DISV = 2 !< DISU6 discretization
    ENUMERATOR :: DISU = 3 !< DISV6 discretization
    ! -- 1D discretizations
    ENUMERATOR :: DIS1D = 101 !< DIS1D6 discretization
    ENUMERATOR :: DISV1D = 102 !< DISV1D6 discretization
    ENUMERATOR :: DISU1D = 103 !< DISU1D6 discretization
    ! -- 2D discretizations
    ENUMERATOR :: DIS2D = 201 !< DIS2D6 discretization
    ENUMERATOR :: DISV2D = 202 !< DISV2D6 discretization
    ENUMERATOR :: DISU2D = 203 !< DISU2D6 discretization
  END ENUM

  !> @brief enumerator used with table objects
  !<
  ENUM, BIND(C)
    ENUMERATOR :: TABLEFT = 0 !< left justified table column
    ENUMERATOR :: TABCENTER = 1 !< centered table column
    ENUMERATOR :: TABRIGHT = 2 !< right justified table column
  END ENUM

  !> @brief enumerator used to define table column data type
  !<
  ENUM, BIND(C)
    ENUMERATOR :: TABSTRING = 0 !< string table data
    ENUMERATOR :: TABUCSTRING = 1 !< upper case string table data
    ENUMERATOR :: TABINTEGER = 2 !< integer table data
    ENUMERATOR :: TABREAL = 3 !< real table data
  END ENUM

  !> @brief enumerator used to define output option
  !<
  ENUM, BIND(C)
    ENUMERATOR :: VSUMMARY = 0 !< write summary output
    ENUMERATOR :: VALL = 1 !< write all simulation notes and warnings
    ENUMERATOR :: VDEBUG = 2 !< write debug output
  END ENUM

  !> @brief enumerator that defines the operating system
  !<
  ENUM, BIND(C)
    ENUMERATOR :: OSUNDEF = 0 !< unknown operating system
    ENUMERATOR :: OSLINUX = 1 !< Linux operating system
    ENUMERATOR :: OSMAC = 2 !< MacOS operating system
    ENUMERATOR :: OSWIN = 3 !< Windows operating system
  END ENUM

  !> @brief enumerator that defines the simulation mode
  !<
  ENUM, BIND(C)
    ENUMERATOR :: MVALIDATE = 0 !< validation mode - do not run time steps
    ENUMERATOR :: MNORMAL = 1 !< normal output mode
    ENUMERATOR :: MRUN = 2 !< run output mode
  END ENUM

  !> @brief enumerator that defines the compiler
  !<
  ENUM, BIND(C)
    ENUMERATOR :: CUNKNOWN = 0 !< unknown compiler
    ENUMERATOR :: CGFORTRAN = 1 !< gfortran compiler
    ENUMERATOR :: CINTEL = 3 !< intel ifort compiler
    ENUMERATOR :: CCRAYFTN = 3 !< cray fortran compiler
  END ENUM

  !> @brief enumerator that defines the cell connection type
  !<
  ENUM, BIND(C)
    ENUMERATOR :: C3D_VERTICAL = 0 !< vertical connection
    ENUMERATOR :: C3D_HORIZONTAL = 1 !< horizontal connection
    ENUMERATOR :: C3D_STAGGERED = 2 !< horizontal connection for a vertically staggered grid
  END ENUM

end module ConstantsModule
