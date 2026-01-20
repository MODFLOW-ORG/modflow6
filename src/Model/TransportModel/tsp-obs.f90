module TspObsModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, MAXOBSTYPES, LENVARNAME
  use BaseDisModule, only: DisBaseType
  use TspIcModule, only: TspIcType
  use ObserveModule, only: ObserveType
  use ObsModule, only: ObsType
  use SimModule, only: count_errors, store_error, &
                       store_error_unit

  implicit none

  private
  public :: TspObsType, tsp_obs_cr

  type, extends(ObsType) :: TspObsType
    ! -- Private members
    type(TspIcType), pointer, private :: ic => null() ! initial conditions
    real(DP), dimension(:), pointer, contiguous, private :: x => null() ! concentration or temperature
    real(DP), dimension(:), pointer, contiguous, private :: flowja => null() ! intercell flows
    character(len=LENVARNAME) :: depvartype = '' !< "concentration" or "temperature"
  contains
    ! -- Public procedures
    procedure, public :: tsp_obs_ar
    procedure, public :: obs_bd => tsp_obs_bd
    procedure, public :: obs_df => tsp_obs_df
    procedure, public :: obs_rp => tsp_obs_rp
    procedure, public :: obs_da => tsp_obs_da
    ! -- Private procedures
    procedure, private :: set_pointers
  end type TspObsType

contains

  !> @brief Create a new TspObsType object
  !!
  !! This routine:
  !!   - creates an observation object
  !!   - allocates pointers
  !!   - initializes values
  !<
  subroutine tsp_obs_cr(obs, inobs, dvt)
    ! -- dummy
    type(TspObsType), pointer, intent(out) :: obs
    integer(I4B), pointer, intent(in) :: inobs
    character(len=LENVARNAME), intent(in) :: dvt !< "concentration" or "temperature"
    !
    allocate (obs)
    call obs%allocate_scalars()
    obs%active = .false.
    obs%inputFilename = ''
    obs%inUnitObs => inobs
    obs%depvartype = dvt
  end subroutine tsp_obs_cr

  !> @brief Allocate and read method for package
  !!
  !!  Method to allocate and read static data for the package.
  !<
  subroutine tsp_obs_ar(this, ic, x, flowja)
    ! -- dummy
    class(TspObsType), intent(inout) :: this
    type(TspIcType), pointer, intent(in) :: ic
    real(DP), dimension(:), pointer, contiguous, intent(in) :: x
    real(DP), dimension(:), pointer, contiguous, intent(in) :: flowja
    !
    ! -- Call ar method of parent class
    call this%obs_ar()
    !
    ! -- set pointers
    call this%set_pointers(ic, x, flowja)
  end subroutine tsp_obs_ar

  !> @brief Define observation object
  !<
  subroutine tsp_obs_df(this, iout, pkgname, filtyp, dis)
    ! -- dummy
    class(TspObsType), intent(inout) :: this
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: filtyp
    class(DisBaseType), pointer :: dis
    ! -- local
    integer(I4B) :: indx
    !
    ! -- Call overridden method of parent class
    call this%ObsType%obs_df(iout, pkgname, filtyp, dis)
    !
    ! -- StoreObsType arguments are: (ObserveType, cumulative, indx);
    !    indx is returned.
    !
    ! -- Store obs type and assign procedure pointer for head observation type
    call this%StoreObsType(trim(adjustl(this%depvartype)), .false., indx)
    this%obsData(indx)%ProcessIdPtr => tsp_process_obs_id
    !
    ! -- Store obs type and assign procedure pointer for flow-ja-face observation type
    call this%StoreObsType('flow-ja-face', .true., indx)
    this%obsData(indx)%ProcessIdPtr => tsp_process_intercell_obs_id
  end subroutine tsp_obs_df

  !> @brief Save observations
  !<
  subroutine tsp_obs_bd(this)
    ! -- dummy
    class(TspObsType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, jaindex, nodenumber
    character(len=100) :: msg
    class(ObserveType), pointer :: obsrv => null()
    !
    call this%obs_bd_clear()
    !
    ! -- iterate through all GWT observations
    if (this%npakobs > 0) then
      do i = 1, this%npakobs
        obsrv => this%pakobs(i)%obsrv
        nodenumber = obsrv%NodeNumber
        jaindex = obsrv%JaIndex
        select case (obsrv%ObsTypeId)
        case ('CONCENTRATION', 'TEMPERATURE')
          call this%SaveOneSimval(obsrv, this%x(nodenumber))
        case ('FLOW-JA-FACE')
          call this%SaveOneSimval(obsrv, this%flowja(jaindex))
        case default
          msg = ' Unrecognized observation type: '//trim(obsrv%ObsTypeId)
          call store_error(msg)
          call store_error_unit(this%inUnitObs)
        end select
      end do
    end if
  end subroutine tsp_obs_bd

  !> @brief If transport model observations need checks, add them here
  !<
  subroutine tsp_obs_rp(this)
    ! -- dummy
    class(TspObsType), intent(inout) :: this
    !
    ! Do GWT (or GWE) observations need any checking? If so, add checks here
  end subroutine tsp_obs_rp

  !> Deallocate memory
  !!
  !! Deallocate memory associated with transport model
  !<
  subroutine tsp_obs_da(this)
    ! -- dummy
    class(TspObsType), intent(inout) :: this
    !
    nullify (this%ic)
    nullify (this%x)
    nullify (this%flowja)
    call this%ObsType%obs_da()
  end subroutine tsp_obs_da

  !> @brief Set pointers needed by the transport OBS package
  !<
  subroutine set_pointers(this, ic, x, flowja)
    ! -- dummy
    class(TspObsType), intent(inout) :: this
    type(TspIcType), pointer, intent(in) :: ic
    real(DP), dimension(:), pointer, contiguous, intent(in) :: x
    real(DP), dimension(:), pointer, contiguous, intent(in) :: flowja
    !
    this%ic => ic
    this%x => x
    this%flowja => flowja
  end subroutine set_pointers

  !> @brief Procedure related to Tsp observations (NOT type-bound)
  !!
  !! Process a specific observation ID
  !<
  subroutine tsp_process_obs_id(obsrv, dis, inunitobs, iout)
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: nn1
    integer(I4B) :: icol, istart, istop
    character(len=LINELENGTH) :: ermsg, string
    !
    ! -- Initialize variables
    string = obsrv%IDstring
    icol = 1
    !
    ! Get node number, with option for ID string to be either node
    ! number or lay, row, column (when dis is structured).
    nn1 = dis%noder_from_string(icol, istart, istop, inunitobs, &
                                iout, string, .false.)
    !
    if (nn1 > 0) then
      obsrv%NodeNumber = nn1
    else
      ermsg = 'Error reading data from ID string'
      call store_error(ermsg)
      call store_error_unit(inunitobs)
    end if
  end subroutine tsp_process_obs_id

  !> @brief Procedure related to Tsp observations (NOT type-bound)
  !!
  !! Process an intercell observation requested by the user
  !<
  subroutine tsp_process_intercell_obs_id(obsrv, dis, inunitobs, iout)
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: nn1, nn2
    integer(I4B) :: icol, istart, istop, jaidx
    character(len=LINELENGTH) :: ermsg, string
    ! formats
70  format('Error: No connection exists between cells identified in text: ', a)
    !
    ! -- Initialize variables
    string = obsrv%IDstring
    icol = 1
    !
    ! Get node number, with option for ID string to be either node
    ! number or lay, row, column (when dis is structured).
    nn1 = dis%noder_from_string(icol, istart, istop, inunitobs, &
                                iout, string, .false.)
    !
    if (nn1 > 0) then
      obsrv%NodeNumber = nn1
    else
      ermsg = 'Error reading data from ID string: '//string(istart:istop)
      call store_error(ermsg)
    end if
    !
    ! Get node number, with option for ID string to be either node
    ! number or lay, row, column (when dis is structured).
    nn2 = dis%noder_from_string(icol, istart, istop, inunitobs, &
                                iout, string, .false.)
    if (nn2 > 0) then
      obsrv%NodeNumber2 = nn2
    else
      ermsg = 'Error reading data from ID string: '//string(istart:istop)
      call store_error(ermsg)
    end if
    !
    ! -- store JA index
    jaidx = dis%con%getjaindex(nn1, nn2)
    if (jaidx == 0) then
      write (ermsg, 70) trim(string)
      call store_error(ermsg)
    end if
    obsrv%JaIndex = jaidx
    !
    if (count_errors() > 0) then
      call store_error_unit(inunitobs)
    end if
  end subroutine tsp_process_intercell_obs_id

end module TspObsModule
