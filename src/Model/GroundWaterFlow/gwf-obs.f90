module GwfObsModule

  use KindModule, only: DP, I4B
  use ConstantsModule, only: LINELENGTH, MAXOBSTYPES
  use BaseDisModule, only: DisBaseType
  use GwfIcModule, only: GwfIcType
  use ObserveModule, only: ObserveType
  use ObsModule, only: ObsType
  use SimModule, only: count_errors, store_error, &
                       store_error_unit
  implicit none

  private
  public :: GwfObsType, gwf_obs_cr

  type, extends(ObsType) :: GwfObsType
    ! -- Private members
    type(GwfIcType), pointer, private :: ic => null() ! initial conditions
    real(DP), dimension(:), pointer, contiguous, private :: x => null() ! head
    real(DP), dimension(:), pointer, contiguous, private :: flowja => null() ! intercell flows

  contains

    ! -- Public procedures
    procedure, public :: gwf_obs_ar
    procedure, public :: obs_bd => gwf_obs_bd
    procedure, public :: obs_df => gwf_obs_df
    procedure, public :: obs_rp => gwf_obs_rp
    procedure, public :: obs_da => gwf_obs_da
    ! -- Private procedures
    procedure, private :: set_pointers
  end type GwfObsType

contains

  !> @brief Create a new GwfObsType object
  !!
  !! Create observation object, allocate pointers, initialize values
  !<
  subroutine gwf_obs_cr(obs, inobs)
    ! -- dummy
    type(GwfObsType), pointer, intent(out) :: obs
    integer(I4B), pointer, intent(in) :: inobs
    !
    allocate (obs)
    call obs%allocate_scalars()
    obs%active = .false.
    obs%inputFilename = ''
    obs%inUnitObs => inobs
  end subroutine gwf_obs_cr

  !> @brief Allocate and read
  !<
  subroutine gwf_obs_ar(this, ic, x, flowja)
    ! -- dummy
    class(GwfObsType), intent(inout) :: this
    type(GwfIcType), pointer, intent(in) :: ic
    real(DP), dimension(:), pointer, contiguous, intent(in) :: x
    real(DP), dimension(:), pointer, contiguous, intent(in) :: flowja
    !
    ! Call ar method of parent class
    call this%obs_ar()
    !
    ! set pointers
    call this%set_pointers(ic, x, flowja)
  end subroutine gwf_obs_ar

  !> @brief Define
  !<
  subroutine gwf_obs_df(this, iout, pkgname, filtyp, dis)
    ! -- dummy
    class(GwfObsType), intent(inout) :: this
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: pkgname
    character(len=*), intent(in) :: filtyp
    class(DisBaseType), pointer :: dis
    ! -- local
    integer(I4B) :: indx
    !
    ! Call overridden method of parent class
    call this%ObsType%obs_df(iout, pkgname, filtyp, dis)
    !
    ! -- StoreObsType arguments are: (ObserveType, cumulative, indx);
    !    indx is returned.
    !
    ! -- Store obs type and assign procedure pointer for head observation type
    call this%StoreObsType('head', .false., indx)
    this%obsData(indx)%ProcessIdPtr => gwf_process_head_drawdown_obs_id
    !
    ! -- Store obs type and assign procedure pointer for drawdown observation type
    call this%StoreObsType('drawdown', .false., indx)
    this%obsData(indx)%ProcessIdPtr => gwf_process_head_drawdown_obs_id
    !
    ! -- Store obs type and assign procedure pointer for flow-ja-face observation type
    call this%StoreObsType('flow-ja-face', .true., indx)
    this%obsData(indx)%ProcessIdPtr => gwf_process_intercell_obs_id
  end subroutine gwf_obs_df

  !> @brief Save obs
  !<
  subroutine gwf_obs_bd(this)
    ! -- dummy
    class(GwfObsType), intent(inout) :: this
    ! -- local
    integer(I4B) :: i, jaindex, nodenumber
    real(DP) :: v
    character(len=100) :: msg
    class(ObserveType), pointer :: obsrv => null()
    !
    call this%obs_bd_clear()
    !
    ! -- iterate through all GWF observations
    if (this%npakobs > 0) then
      do i = 1, this%npakobs
        obsrv => this%pakobs(i)%obsrv
        nodenumber = obsrv%NodeNumber
        jaindex = obsrv%JaIndex
        select case (obsrv%ObsTypeId)
        case ('HEAD')
          call this%SaveOneSimval(obsrv, this%x(nodenumber))
        case ('DRAWDOWN')
          v = this%ic%strt(nodenumber) - this%x(nodenumber)
          call this%SaveOneSimval(obsrv, v)
        case ('FLOW-JA-FACE')
          call this%SaveOneSimval(obsrv, this%flowja(jaindex))
        case default
          msg = 'Error: Unrecognized observation type: '//trim(obsrv%ObsTypeId)
          call store_error(msg)
        end select
      end do
      !
      ! -- write summary of error messages
      if (count_errors() > 0) then
        call store_error_unit(this%inUnitObs)
      end if
    end if
  end subroutine gwf_obs_bd

  !> @brief Do GWF observations need any checking? If so, add checks here
  !<
  subroutine gwf_obs_rp(this)
    ! -- dummy
    class(GwfObsType), intent(inout) :: this
    !
    ! Do GWF observations need any checking? If so, add checks here
  end subroutine gwf_obs_rp

  !> @brief Deallocate memory
  !<
  subroutine gwf_obs_da(this)
    ! -- dummy
    class(GwfObsType), intent(inout) :: this
    !
    nullify (this%ic)
    nullify (this%x)
    nullify (this%flowja)
    call this%ObsType%obs_da()
  end subroutine gwf_obs_da

  !> @brief Set pointers
  !<
  subroutine set_pointers(this, ic, x, flowja)
    ! -- dummy
    class(GwfObsType), intent(inout) :: this
    type(GwfIcType), pointer, intent(in) :: ic
    real(DP), dimension(:), pointer, contiguous, intent(in) :: x
    real(DP), dimension(:), pointer, contiguous, intent(in) :: flowja
    !
    this%ic => ic
    this%x => x
    this%flowja => flowja
  end subroutine set_pointers

  ! -- Procedures related to GWF observations (NOT type-bound)

  !> @brief Calculate drawdown observation when requested
  !<
  subroutine gwf_process_head_drawdown_obs_id(obsrv, dis, inunitobs, iout)
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
  end subroutine gwf_process_head_drawdown_obs_id

  !> @brief Process flow between two cells when requested
  !<
  subroutine gwf_process_intercell_obs_id(obsrv, dis, inunitobs, iout)
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
  end subroutine gwf_process_intercell_obs_id

end module GwfObsModule
