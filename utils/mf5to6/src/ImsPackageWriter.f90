module ImsPackageWriterModule

  use ConstantsModule, only: DZERO, MAXCHARLEN
  use ConstantsPHMFModule, only: FCINPUT, LENCTYPE
  use DE4MODULE, only: itmxd4 => ITMX, MUTD4, HCLOSEDE4, IPRD4, ACCLDE4, &
                       ifreqd4 => IFREQ
  use GLOBAL, only: NLAY, NROW, NCOL
  use GlobalVariablesModule, only: ilgr
  use GMGMODULE, only: iitergmg => IITER, HCLOSEGMG, RCLOSEGMG, DAMPGMG, &
                       IADAMPGMG, RELAXGMG, mxitergmg
  use GWFNWTMODULE, only: Ftol, mxiterNwt => ITER1, IPRNWT, Tol, &
                          MxIterInner, HCloseLinear, RCloseLinear, Theta, &
                          Akappa, GammaNwt => Gamma, Amomentum, &
                          Linmeth, IPRNWT, IBOTAV, &
                          Iter1Nwt => iter1, IFDPARAM, iaclNwt, norderNwt, &
                          levelNwt, northNwt, MaxBackIterNwt, IRedSysNwt, &
                          IDropTolNwt, BackFlagNwt, BackTolNwt, &
                          BackReduceNwt, EpsrnNwt, RrctolNwt
  use InputOutputModule, only: GetUnit
  use PackageWriterModule, only: PackageWriterType
  use PCGMODULE, only: mxiterpcg, MUTPCG, HCLOSEPCG, IPRPCG, iter1pcg => ITER1, &
                       RCLOSEPCG, RELAXPCG, DAMPPCG, npcondpcg => NPCOND
  use pcgn, only: pcgndat
  use SimPHMFModule, only: ustop, store_error, store_note, store_warning
  use SIPMODULE, only: mxitersip, hclosesip => HCLOSE, IPRSIP, acclsip => ACCL
  use UtilitiesModule, only: GreaterOf, GreatestOf
  !use XMDMODULE, only: LevelXmd => LEVEL, EpsrnXmd => EPSRN, NorthXmd => NORTH

  integer, parameter :: MODERATETHRESHOLD = 100000
  integer, parameter :: COMPLEXTHRESHOLD = 500000

  type, extends(PackageWriterType) :: ImsPackageWriterType
    ! MF-2005 solver
    character(len=4) :: mf2005_solver = ''
    ! IMS input variables and default values
    ! Options
    character(len=7) :: print_option = 'SUMMARY'
    character(len=8) :: complexity = 'SIMPLE'
    ! Nonlinear solver
    real             :: outer_dvclose = 0.01
    integer          :: outer_maximum = 100
    character(len=6) :: under_relaxation = 'NONE'
    double precision :: UNDER_RELAXATION_THETA = DZERO
    double precision :: UNDER_RELAXATION_KAPPA = DZERO
    double precision :: UNDER_RELAXATION_GAMMA = DZERO
    double precision :: UNDER_RELAXATION_MOMENTUM = DZERO
    character(len=4) :: linear_solver = 'PCGU'
!    character(len=4) :: linear_solver = 'XMD'
    ! Linear solvers
    real             :: inner_dvclose = 0.001
    real             :: inner_rclose = 0.1
    double precision :: preconditioner_drop_tolerance = -0.001d0
    integer          :: inner_maximum = 100
    integer          :: preconditioner_levels = -7
    integer          :: number_orthogonalizations = 2
    character(len=8) :: linear_acceleration = 'CG'
!    character(len=8) :: linear_acceleration = 'BICGSTAB'
    character(len=4) :: REORDERING_METHOD = 'NONE'
    double precision :: RELAXATION_FACTOR = DZERO
    character(len=6) :: preconditioner = 'JACOBI'
!    character(len=6) :: preconditioner = 'NONE'
    integer          :: backtracking_number = 0
    double precision :: backtracking_tolerance = 1.0d4
    double precision :: backtracking_reduction_factor = 0.2d0
    double precision :: backtracking_residual_limit = 100.0d0
  contains
    procedure, public :: MyType
    procedure, public :: ProcessAllocate
    procedure, public :: ProcessStressLoop
    procedure, public :: WriteFile
    procedure, private :: DefineVariables
  end type ImsPackageWriterType

  integer, parameter :: MINOUTER = 20

contains

  subroutine DefineVariables(this, igrid)
    implicit none
    ! dummy
    class(ImsPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    character(len=MAXCHARLEN) :: msg
    integer :: nnodes
    !
    ! Set default complexity, linear solver, linear acceleration
    ! if it seems it would help.
    nnodes = NLAY * NROW * NCOL
    if (nnodes > COMPLEXTHRESHOLD) then
      this%complexity = 'COMPLEX'
      this%linear_acceleration = 'BICGSTAB'
    elseif (nnodes > MODERATETHRESHOLD) then
      this%complexity = 'MODERATE'
      this%linear_acceleration = 'CG'
    else
      ! Defaults
      this%complexity = 'SIMPLE'
      this%linear_solver = 'PCGU'
      this%linear_acceleration = 'CG'
    endif
    if (ilgr > 0) then
      this%linear_acceleration = 'BICGSTAB'
      this%RELAXATION_FACTOR = 0.97D0
    endif
    !
    select case (this%mf2005_solver)
    case ('DE4')
      this%outer_dvclose = HCLOSEDE4
      if (ifreqd4==1) then
        this%outer_maximum = GreaterOf(this%outer_maximum, MINOUTER)
      else
        this%outer_maximum = GreatestOf(this%outer_maximum, MINOUTER, itmxd4)
        if (ACCLDE4 < 1.0 .and. this%linear_solver=='PCGU') then
          msg = 'In DE4 package input, ACCL < 1. You may need to ' // &
                'use under relaxation in IMS input.'
          call store_note(msg)
        endif
      endif
    case ('GMG')
      this%outer_dvclose = HCLOSEGMG
      this%outer_maximum = GreatestOf(this%outer_maximum, MINOUTER, mxitergmg)
      this%inner_rclose = RCLOSEGMG
      this%inner_maximum = iitergmg
      if(ilgr == 0) this%RELAXATION_FACTOR = RELAXGMG
      if (mxitergmg > 1 .and. this%linear_solver=='PCGU') then
        if (DAMPGMG < 1.0) then
          msg = 'In GMG package input, DAMP < 1. You may need to ' // &
                'use under relaxation in IMS input.'
          call store_note(msg)
        elseif (IADAMPGMG /= 0) then
          msg = 'In GMG package input, IADAMP is not 0. You may need to ' // &
                'use under relaxation in IMS input.'
          call store_note(msg)
        endif
      endif
    case ('PCG')
      this%outer_dvclose = HCLOSEPCG
      this%outer_maximum = GreatestOf(this%outer_maximum, MINOUTER, mxiterpcg)
      this%inner_maximum = iter1pcg
      this%inner_rclose = RCLOSEPCG
      if(ilgr == 0) this%RELAXATION_FACTOR = RELAXPCG
      this%preconditioner_drop_tolerance = -1.0d0
      this%preconditioner_levels = -1
      if (mxiterpcg > 1 .and. this%linear_solver=='PCGU') then
        if (DAMPPCG < 1.0 .or. (npcondpcg==1 .and. RELAXPCG < 1.0)) then
          msg = 'In PCG package input, damping and/or relaxation is used. You may ' // &
                'need to use under relaxation in IMS input.'
          call store_note(msg)
        endif
      endif
    case ('PCGN')
      this%outer_dvclose = PCGNDAT(igrid)%HCLOSE
      this%outer_maximum = GreatestOf(this%outer_maximum, MINOUTER, PCGNDAT(igrid)%MO_ITER)
      this%inner_maximum = PCGNDAT(igrid)%MI_ITER
      this%inner_rclose = PCGNDAT(igrid)%RCLOSE
      if(ilgr == 0) this%RELAXATION_FACTOR = PCGNDAT(igrid)%DRELAX
      this%preconditioner_levels = PCGNDAT(igrid)%FILL
      if (PCGNDAT(igrid)%MO_ITER > 1 .and. this%linear_solver=='PCGU') then
        if (PCGNDAT(igrid)%SAV_DAMP < 1.0) then
          msg = 'In PCGN package input, DAMP < 1. You may need to ' // &
                'use under relaxation in IMS input.'
          call store_note(msg)
        elseif (PCGNDAT(igrid)%DAMP_A /= 0) then
          msg = 'In PCGN package input, ADAMP is not 0. You may need to ' // &
                'use under relaxation in IMS input.'
          call store_note(msg)
        endif
      endif
    case ('SIP')
      this%outer_dvclose = hclosesip
      this%outer_maximum = GreatestOf(this%outer_maximum, MINOUTER, mxitersip)
      if (acclsip < 1.0 .and. this%linear_solver=='PCGU') then
        msg = 'In SIP package input, ACCL < 1. You may need to ' // &
              'use under relaxation in IMS input.'
        call store_note(msg)
      endif
    case ('NWT')
      ! Nonlinear parameters
      this%outer_dvclose = Tol
      this%outer_maximum = GreatestOf(this%outer_maximum, MINOUTER, mxiterNwt)
      !this%RELAXATION_FACTOR = RelaxNwt
      this%under_relaxation = 'DBD'
      this%UNDER_RELAXATION_THETA = Theta
      this%UNDER_RELAXATION_KAPPA = Akappa
      this%UNDER_RELAXATION_GAMMA = GammaNwt
      this%UNDER_RELAXATION_MOMENTUM = Amomentum
      if (BackFlagNwt > 0) then
!       Code from SUBROUTINE GWF2NWT1AR:
!       if (Btrack > 0) then
!         BackFlagNwt = Btrack
!         MaxBackIterNwt = Numtrack
!         BackTolNwt = Btol
!         BackReduceNwt = Breduc
!       endif
        this%backtracking_number = MaxBackIterNwt
        this%backtracking_tolerance = BackTolNwt
        this%backtracking_reduction_factor = BackReduceNwt
      endif
      ! Linear parameters
      this%inner_maximum = MxIterInner
      if (RrctolNwt > 0.0d0) then
        this%inner_rclose = RrctolNwt
      endif
      select case (iaclNwt)
      case (0)
        this%linear_acceleration = 'BICGSTAB'
      case (1)
        this%linear_acceleration = 'BICGSTAB'
      case (2)
        this%linear_acceleration = 'BICGSTAB'
      end select
      select case (norderNwt)
      case (0)
        this%reordering_method = 'NONE'
      case (1)
        this%reordering_method = 'RCM'
      case (2)
        this%reordering_method = 'MD'
      end select
      this%preconditioner_levels = levelNwt
      if (IDropTolNwt == 0) then
        this%preconditioner_drop_tolerance = 0.0d0
      else
        this%preconditioner_drop_tolerance = EpsrnNwt
      endif
      this%number_orthogonalizations = northNwt
    case ('DFLT')
      msg = 'Solver input file not defined. ' // &
          ' Default input for IMS solver will be provided.'
      call store_warning(msg)
    case default
      msg = 'Solver not recognized: ' // trim(this%mf2005_solver)
      call store_warning(msg)
    end select
    !
    ! Adjust inner_dvclose if it's too large
    if (associated(HCloseLinear)) then
      this%inner_dvclose = HCloseLinear
    endif
    if (this%linear_solver == 'XMD') then
      if (this%inner_dvclose > (this%outer_dvclose * 0.1)) then
        this%inner_dvclose = 0.1 * this%outer_dvclose
      endif
    elseif (this%linear_solver == 'PCGU') then
      if (this%inner_dvclose > (this%outer_dvclose * 0.01)) then
        this%inner_dvclose = 0.1 * this%outer_dvclose
      endif
    endif
    !!
    !! Set default complexity, linear solver, linear acceleration
    !! if it seems it would help.
    !nnodes = NLAY * NROW * NCOL
    !if (nnodes > COMPLEXTHRESHOLD) then
    !  this%complexity = 'COMPLEX'
    !  !this%linear_solver = 'XMD'
    !  this%linear_acceleration = 'BICGSTAB'
    !elseif (nnodes > MODERATETHRESHOLD) then
    !  this%complexity = 'MODERATE'
    !  !this%linear_solver = 'XMD'
    !  this%linear_acceleration = 'CG'
    !else
    !  ! Defaults
    !  this%complexity = 'SIMPLE'
    !  this%linear_solver = 'PCGU'
    !  this%linear_acceleration = 'CG'
    !endif
    !
    return
  end subroutine DefineVariables

  subroutine WriteFile(this, igrid)
    ! Write the IMS file
    implicit none
    ! dummy
    class(ImsPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! local
    integer :: iu
    character(len=MAXCHARLEN) :: fname
    character(len=4) :: ftype = 'IMS6'
    ! formats
    1 format()
    10 format(a)
    20 format(2x,a,2x,i0)
    30 format(2x,a,2x,a)
    40 format(2x,a,2x,g14.7)
    !
    ! Initialize IMS file
    this%Active = .true.
    this%fileobj%FCode = FCINPUT
    fname = trim(this%ModelBasename) // '.ims'
    call this%InitializeFile(fname, ftype, this%PackageName)
    !
    ! Define IMS input variables
    call this%DefineVariables(igrid)
    !
    !
    iu = this%fileobj%IUnit
    write(iu,1)
    write(iu,10)'BEGIN Options'
    write(iu,30)'PRINT_OPTION', this%print_option
    write(iu,30)'COMPLEXITY', trim(this%complexity)
    write(iu,10)'END Options'
    write(iu,1)
    write(iu,10)'BEGIN Nonlinear'
    write(iu,40)'OUTER_DVCLOSE', this%outer_dvclose
    write(iu,20)'OUTER_MAXIMUM', this%outer_maximum
    write(iu,30)'UNDER_RELAXATION', trim(this%under_relaxation)
!    if (this%linear_solver /= 'PCGU') then
!      write(iu,30)'LINEAR_SOLVER',trim(this%linear_solver)
!    endif
    write(iu,40)'UNDER_RELAXATION_THETA', this%UNDER_RELAXATION_THETA
    write(iu,40)'UNDER_RELAXATION_KAPPA', this%UNDER_RELAXATION_KAPPA
    write(iu,40)'UNDER_RELAXATION_GAMMA', this%UNDER_RELAXATION_GAMMA
    write(iu,40)'UNDER_RELAXATION_MOMENTUM', this%UNDER_RELAXATION_MOMENTUM
    if (this%backtracking_number > 0) then
      write(iu,20)'BACKTRACKING_NUMBER', this%backtracking_number
      write(iu,40)'BACKTRACKING_TOLERANCE', this%backtracking_tolerance
      write(iu,40)'BACKTRACKING_REDUCTION_FACTOR', this%backtracking_reduction_factor 
      write(iu,40)'BACKTRACKING_RESIDUAL_LIMIT', this%backtracking_residual_limit
    endif
    write(iu,10)'END Nonlinear'
    write(iu,1)
    write(iu,10)'BEGIN LINEAR'
    write(iu,20)'INNER_MAXIMUM', this%inner_maximum
    write(iu,40)'INNER_DVCLOSE', this%inner_dvclose
    write(iu,40)'INNER_RCLOSE', this%inner_rclose
    write(iu,30)'LINEAR_ACCELERATION', trim(this%linear_acceleration)
    if (this%RELAXATION_FACTOR > DZERO) then
      write(iu,40)'RELAXATION_FACTOR', this%RELAXATION_FACTOR
    endif
    if (this%preconditioner_levels > 0) then
      write(iu,20)'PRECONDITIONER_LEVELS', this%preconditioner_levels
    endif
    if (this%preconditioner_drop_tolerance > DZERO) then
      write(iu,40)'PRECONDITIONER_DROP_TOLERANCE', this%preconditioner_drop_tolerance
    endif
    write(iu,20)'NUMBER_ORTHOGONALIZATIONS', this%number_orthogonalizations
!    write(iu,40)'SCALING_METHOD',this  !  <-- NEED variable
    write(iu,30)'REORDERING_METHOD', this%REORDERING_METHOD
    write(iu,10)'END LINEAR'
    call this%CloseFile()
    !
    return
  end subroutine WriteFile

  subroutine ProcessAllocate(this, igrid)
    class(ImsPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! nothing to allocate
    return
  end subroutine ProcessAllocate

  subroutine ProcessStressLoop(this, igrid)
    class(ImsPackageWriterType) :: this
    integer, intent(in) :: igrid
    ! no stresses to process
    return
  end subroutine ProcessStressLoop

  function MyType(this) result (ctype)
    ! dummy 
    class(ImsPackageWriterType) :: this
    character(len=LENCTYPE) :: ctype
    !
    ctype = 'ImsPackageWriterType'
    !
    return
  end function MyType
  
end module ImsPackageWriterModule
