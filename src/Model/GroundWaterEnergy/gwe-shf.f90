!> @brief This module contains methods for calculating sensible heat flux
!!
!! This module contains the methods used to calculate the sensible heat flux
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

module SensHeatModule
  use ConstantsModule, only: LINELENGTH, LENMEMPATH, DZERO, LENVARNAME
  use KindModule, only: I4B, DP
  use MemoryManagerModule, only: mem_setptr
  use MemoryHelperModule, only: create_mem_path
  use TimeSeriesLinkModule, only: TimeSeriesLinkType
  use TimeSeriesManagerModule, only: TimeSeriesManagerType, tsmanager_cr
  use SimModule, only: store_error
  use SimVariablesModule, only: errmsg
  use PbstBaseModule, only: PbstBaseType, pbstbase_da

  implicit none

  private

  public :: ShfType
  public :: shf_cr

  character(len=16) :: text = '          SHF'

  type, extends(PbstBaseType) :: ShfType

    real(DP), pointer :: rhoa => null() !< ABC density of air
    real(DP), pointer :: cpa => null() !< ABC heat capacity of air
    real(DP), pointer :: cd => null() !< ABC drag coefficient
    real(DP), dimension(:), pointer, contiguous :: wspd => null() !< ABC wind speed
    real(DP), dimension(:), pointer, contiguous :: tatm => null() !< ABC temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: patm => null() !< ABC atmospheric pressure
    real(DP), dimension(:), pointer, contiguous :: ea => null() !< ABC temperature of the atmosphere
    real(DP), dimension(:), pointer, contiguous :: ew => null() !< ABC atmospheric pressure

  contains

    procedure :: da => shf_da
    procedure :: pbst_ar => shf_ar_set_pointers
    procedure, public :: shf_cq

  end type ShfType

contains

  !> @brief Create a new ShfType object
  !!
  !! Create a new sensible heat flux (ShfType) object. Initially for use with
  !! the SFE-ABC package.
  !<
  subroutine shf_cr(this, name_model, inunit, iout, ncv)
    ! -- dummy
    type(ShfType), pointer, intent(out) :: this
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    allocate (this)
    call this%init(name_model, 'SHF', 'SHF', inunit, iout, ncv)
    this%text = text
  end subroutine shf_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the ABC
  !! package for access by SHF.
  !<
  subroutine shf_ar_set_pointers(this)
    ! -- dummy
    class(ShfType) :: this
    ! -- local
    character(len=LENMEMPATH) :: abcMemoryPath
    !
    ! -- print a message noting that the SHF utility is active
    write (this%iout, '(a)') &
      'SHF -- SENSIBLE HEAT WILL BE INCLUDED IN THE ATMOSPHERIC BOUNDARY '// &
      'CONDITIONS FOR THE STREAMFLOW ENERGY TRANSPORT PACKAGE'
    !
    ! -- set pointers to variables hosted in the ABC package
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%rhoa, 'RHOA', abcMemoryPath)
    call mem_setptr(this%cpa, 'CPA', abcMemoryPath)
    call mem_setptr(this%cd, 'CD', abcMemoryPath)
    call mem_setptr(this%wspd, 'WSPD', abcMemoryPath)
    call mem_setptr(this%tatm, 'TATM', abcMemoryPath)
    call mem_setptr(this%patm, 'PATM', abcMemoryPath)
    call mem_setptr(this%ea, 'EA', abcMemoryPath)
    call mem_setptr(this%ew, 'EW', abcMemoryPath)
    !
    ! -- create time series manager
    call tsmanager_cr(this%tsmanager, this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
  end subroutine shf_ar_set_pointers

  !> @brief Calculate Sensible Heat Flux
  !!
  !! Calculate and return the sensible heat flux for one reach
  !<
  subroutine shf_cq(this, ifno, tstrm, shflx, lhflx)
    ! -- dummy
    class(ShfType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(in) :: tstrm !< temperature of the stream reach
    real(DP), intent(inout) :: shflx !< calculated sensible heat flux amount
    real(DP), optional, intent(in) :: lhflx !< latent heat flux
    ! -- local
    real(DP) :: shf_const
    real(DP) :: br
    !
    ! -- calculate sensible heat flux using HGS equation
    if (present(lhflx)) then
      br = 0.00061_DP * this%patm(ifno) * &
           ((tstrm - this%tatm(ifno)) / (this%ew(ifno) - this%ea(ifno)))
      shflx = br * lhflx
    else
      shf_const = this%cd * this%cpa * this%rhoa
      shflx = shf_const * this%wspd(ifno) * (this%tatm(ifno) - tstrm)
    end if
  end subroutine shf_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine shf_da(this)
    ! -- dummy
    class(ShfType) :: this
    !
    ! -- nullify pointers to other package variables
    nullify (this%rhoa)
    nullify (this%cpa)
    nullify (this%cd)
    nullify (this%wspd)
    nullify (this%tatm)
    nullify (this%patm)
    nullify (this%ea)
    nullify (this%ew)
    !
    ! -- deallocate parent
    call pbstbase_da(this)
  end subroutine shf_da

end module SensHeatModule
