!> @brief This module contains methods for calculating the shortwave radiation heat flux
!!
!! This module contains the methods used to calculate the shortwave radiation heat flux
!! for surface-water boundaries, like streams and lakes.  In its current form,
!! this class acts like a package to a package, similar to the TVK package that
!! can be invoked from the NPF package.  Once this package is completed in its
!! prototyped form, it will likely be moved around.
!<

module ShortwaveModule
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

  public :: SwrType
  public :: swr_cr

  character(len=16) :: text = '          SWR'

  type, extends(PbstBaseType) :: SwrType

    real(DP), dimension(:), pointer, contiguous :: solr => null() !< solar radiation
    real(DP), dimension(:), pointer, contiguous :: shd => null() !< shade fraction
    real(DP), dimension(:), pointer, contiguous :: swrefl => null() !< shortwave reflectance of water surface

  contains

    procedure :: da => swr_da
    procedure :: pbst_ar => swr_ar_set_pointers
    procedure, public :: swr_cq

  end type SwrType

contains

  !> @brief Create a new SwrType object
  !!
  !! Create a new shortwave radiation flux (SwrType) object. Initially for use with
  !! the SFE package.
  !<
  subroutine swr_cr(this, name_model, inunit, iout, ncv)
    ! -- dummy
    type(SwrType), pointer, intent(out) :: this
    character(len=*), intent(in) :: name_model
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    integer(I4B), target, intent(in) :: ncv
    !
    allocate (this)
    call this%init(name_model, 'SWR', 'SWR', inunit, iout, ncv)
    this%text = text
  end subroutine swr_cr

  !> @brief Announce package and set pointers to variables
  !!
  !! Announce package version and set array and variable pointers from the ABC
  !! package for access by SWR.
  !<
  subroutine swr_ar_set_pointers(this)
    ! -- dummy
    class(SwrType) :: this
    ! -- local
    character(len=LENMEMPATH) :: abcMemoryPath
    !
    ! -- print a message noting that the SWR utility is active
    write (this%iout, '(a)') &
      'SWR -- SHORTWAVE RADIATION WILL BE INCLUDED IN THE ATMOSPHERIC '// &
      'BOUNDARY CONDITIONS FOR THE STREAMFLOW ENERGY TRANSPORT PACKAGE'
    !
    ! -- Set pointers to variables hosted in the ABC package
    abcMemoryPath = create_mem_path(this%name_model, 'ABC')
    call mem_setptr(this%solr, 'SOLR', abcMemoryPath)
    call mem_setptr(this%shd, 'SHD', abcMemoryPath)
    call mem_setptr(this%swrefl, 'SWREFL', abcMemoryPath)
    !
    ! -- create time series manager
    call tsmanager_cr(this%tsmanager, this%iout, &
                      removeTsLinksOnCompletion=.true., &
                      extendTsToEndOfSimulation=.true.)
  end subroutine swr_ar_set_pointers

  !> @brief Calculate Shortwave Radiation Heat Flux
  !!
  !! Calculate and return the shortwave radiation heat flux for one reach
  !<
  subroutine swr_cq(this, ifno, swrflx)
    ! -- dummy
    class(SwrType), intent(inout) :: this
    integer(I4B), intent(in) :: ifno !< stream reach integer id
    real(DP), intent(inout) :: swrflx !< calculated shortwave radiation heat flux amount
    !
    ! -- calculate shortwave radiation heat flux (version: user input of sol rad data)
    swrflx = (1 - this%shd(ifno)) * (1 - this%swrefl(ifno)) * this%solr(ifno)
  end subroutine swr_cq

  !> @brief Deallocate package memory
  !!
  !! Deallocate TVK package scalars and arrays.
  !<
  subroutine swr_da(this)
    ! -- dummy
    class(SwrType) :: this
    !
    ! -- deallocate time series
    nullify (this%shd)
    nullify (this%swrefl)
    nullify (this%solr)
    !
    ! -- deallocate parent
    call pbstbase_da(this)
  end subroutine swr_da

end module ShortwaveModule
