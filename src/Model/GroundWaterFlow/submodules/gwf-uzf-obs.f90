!> @brief UZF package observations
!<
submodule(UzfModule) UzfModuleObs
contains

  !> @brief Return true because uzf package supports observations
  !!
  !! Overrides BndType%bnd_obs_supported
  !<
  module procedure uzf_obs_supported
  uzf_obs_supported = .true.
  end procedure uzf_obs_supported

  !> @brief Implements bnd_df_obs
  !!
  !! Store observation type supported by uzf package.
  !! Overrides BndType%bnd_df_obs
  !<
  module procedure uzf_df_obs
  ! -- local
  integer(I4B) :: indx
  !
  call this%obs%StoreObsType('uzf-gwrch', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for discharge observation type.
  call this%obs%StoreObsType('uzf-gwd', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for discharge observation type.
  call this%obs%StoreObsType('uzf-gwd-to-mvr', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for gwet observation type.
  call this%obs%StoreObsType('uzf-gwet', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for infiltration observation type.
  call this%obs%StoreObsType('infiltration', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for from mover observation type.
  call this%obs%StoreObsType('from-mvr', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for rejected infiltration observation type.
  call this%obs%StoreObsType('rej-inf', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for rejected infiltration to mover observation type.
  call this%obs%StoreObsType('rej-inf-to-mvr', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for uzet observation type.
  call this%obs%StoreObsType('uzet', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for storage observation type.
  call this%obs%StoreObsType('storage', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for net infiltration observation type.
  call this%obs%StoreObsType('net-infiltration', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  !
  !    for water-content observation type.
  call this%obs%StoreObsType('water-content', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => uzf_process_obsID
  end procedure uzf_df_obs

  !> @brief Process each observation
  !!
  !! Only done the first stress period since boundaries are fixed for the
  !! simulation
  !<
  module procedure uzf_rp_obs
  ! -- modules
  use TdisModule, only: kper
  ! -- local
  integer(I4B) :: i
  integer(I4B) :: j
  integer(I4B) :: n
  integer(I4B) :: nn
  integer(I4B) :: iuzid
  real(DP) :: obsdepth
  real(DP) :: dmax
  character(len=LENBOUNDNAME) :: bname
  class(ObserveType), pointer :: obsrv => null()
  !
60 format('Invalid node number in OBS input: ', i0)
  !
  if (kper == 1) then
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      !
      ! -- get node number 1
      nn = obsrv%NodeNumber
      if (nn == NAMEDBOUNDFLAG) then
        bname = obsrv%FeatureName
        !
        ! -- Observation location(s) is(are) based on a boundary name.
        !    Iterate through all boundaries to identify and store
        !    corresponding index(indices) in bound array.
        do n = 1, this%nodes
          if (this%boundname(n) == bname) then
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(n)
            if (obsrv%indxbnds_count == 1) then
              !
              ! -- Define intPak1 so that obs_theta is stored (for first uzf
              !    cell if multiple cells share the same boundname).
              obsrv%intPak1 = n
            end if
          end if
        end do
      else
        !
        ! -- get node number
        nn = obsrv%NodeNumber
        !
        ! -- put nn (a value meaningful only to UZF) in intPak1
        obsrv%intPak1 = nn
        ! -- check that node number is valid; call store_error if not
        if (nn < 1 .or. nn > this%nodes) then
          write (errmsg, 60) nn
          call store_error(errmsg)
        else
          obsrv%BndFound = .true.
        end if
        obsrv%CurrentTimeStepEndValue = DZERO
        call obsrv%AddObsIndex(nn)
      end if
      !
      ! -- catch non-cumulative observation assigned to observation defined
      !    by a boundname that is assigned to more than one element
      if (obsrv%ObsTypeId == 'WATER-CONTENT') then
        n = obsrv%indxbnds_count
        if (n /= 1) then
          write (errmsg, '(a,3(1x,a))') &
            trim(adjustl(obsrv%ObsTypeId)), 'for observation', &
            trim(adjustl(obsrv%Name)), &
            'must be assigned to a UZF cell with a unique boundname.'
          call store_error(errmsg, terminate=.TRUE.)
        end if
        !
        ! -- check WATER-CONTENT depth
        obsdepth = obsrv%Obsdepth
        !
        ! -- put obsdepth (a value meaningful only to UZF) in dblPak1
        obsrv%dblPak1 = obsdepth
        !
        ! -- determine maximum cell depth
        ! -- This is presently complicated for landflag = 1 cells and surfdep
        !    greater than zero.  In this case, celtop is dis%top - surfdep.
        iuzid = obsrv%intPak1
        dmax = this%uzfobj%celtop(iuzid) - this%uzfobj%celbot(iuzid)
        ! -- check that obs depth is valid; call store_error if not
        ! -- need to think about a way to put bounds on this depth
        ! -- Also, an observation depth of 0.0, whether a landflag == 1 object
        ! -- or a subsurface object, is not legit since this would be at a
        ! -- a layer interface and therefore a discontinuity.
        if (obsdepth <= DZERO .or. obsdepth > dmax) then
          write (errmsg, '(a,3(1x,a),1x,g0,1x,a,1x,g0,a)') &
            trim(adjustl(obsrv%ObsTypeId)), 'for observation', &
            trim(adjustl(obsrv%Name)), 'specified depth (', obsdepth, &
            ') must be greater than 0.0 and less than ', dmax, '.'
          call store_error(errmsg)
        end if
      else
        do j = 1, obsrv%indxbnds_count
          nn = obsrv%indxbnds(j)
          if (nn < 1 .or. nn > this%maxbound) then
            write (errmsg, '(a,2(1x,a),1x,i0,1x,a,1x,i0,a)') &
              trim(adjustl(obsrv%ObsTypeId)), 'uzfno must be greater than 0 ', &
              'and less than or equal to', this%maxbound, &
              '(specified value is ', nn, ').'
            call store_error(errmsg)
          end if
        end do
      end if
    end do
    !
    ! -- evaluate if there are any observation errors
    if (count_errors() > 0) then
      call store_error_unit(this%inunit)
    end if
  end if
  end procedure uzf_rp_obs

  !> @brief Calculate observations this time step and call ObsType%SaveOneSimval
  !! for each UzfType observation
  !<
  module procedure uzf_bd_obs
  ! -- local
  integer(I4B) :: i
  integer(I4B) :: ii
  integer(I4B) :: n
  real(DP) :: v
  type(ObserveType), pointer :: obsrv => null()
  !
  call this%uzf_solve(reset_state=.false.)
  !
  ! Write simulated values for all uzf observations
  if (this%obs%npakobs > 0) then
    call this%obs%obs_bd_clear()
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      do ii = 1, obsrv%indxbnds_count
        n = obsrv%indxbnds(ii)
        v = DNODATA
        select case (obsrv%ObsTypeId)
        case ('UZF-GWRCH')
          v = this%rch(n)
        case ('UZF-GWD')
          v = this%gwd(n)
          if (v > DZERO) then
            v = -v
          end if
        case ('UZF-GWD-TO-MVR')
          if (this%imover == 1) then
            v = this%gwdtomvr(n)
            if (v > DZERO) then
              v = -v
            end if
          end if
        case ('UZF-GWET')
          if (this%igwetflag > 0) then
            v = this%gwet_pvar(n)
            if (v > DZERO) then
              v = -v
            end if
          end if
        case ('INFILTRATION')
          v = this%appliedinf(n)
        case ('FROM-MVR')
          if (this%imover == 1) then
            v = this%pakmvrobj%get_qfrommvr(n)
          end if
        case ('REJ-INF')
          v = this%rejinf(n)
          if (v > DZERO) then
            v = -v
          end if
        case ('REJ-INF-TO-MVR')
          if (this%imover == 1) then
            v = this%rejinftomvr(n)
            if (v > DZERO) then
              v = -v
            end if
          end if
        case ('UZET')
          if (this%ietflag /= 0) then
            v = this%uzet(n)
            if (v > DZERO) then
              v = -v
            end if
          end if
        case ('STORAGE')
          v = -this%qsto(n)
        case ('NET-INFILTRATION')
          v = this%infiltration(n)
        case ('WATER-CONTENT')
          v = this%uzfobj%get_water_content_at_depth(n, obsrv%obsDepth)
        case default
          errmsg = 'Unrecognized observation type: '//trim(obsrv%ObsTypeId)
          call store_error(errmsg)
        end select
        call this%obs%SaveOneSimval(obsrv, v)
      end do
    end do
    !
    ! -- write summary of error messages
    if (count_errors() > 0) then
      call this%parser%StoreErrorUnit()
    end if
  end if
  end procedure uzf_bd_obs

  !> @brief This procedure is pointed to by ObsDataType%ProcesssIdPtr
  !!
  !! Process the ID string of an observation definition for UZF-package
  !! observations
  !<
  subroutine uzf_process_obsID(obsrv, dis, inunitobs, iout)
    ! -- .
    ! -- dummy
    type(ObserveType), intent(inout) :: obsrv
    class(DisBaseType), intent(in) :: dis
    integer(I4B), intent(in) :: inunitobs
    integer(I4B), intent(in) :: iout
    ! -- local
    integer(I4B) :: n, nn
    real(DP) :: obsdepth
    integer(I4B) :: icol, istart, istop, istat
    real(DP) :: r
    character(len=LINELENGTH) :: string
    ! formats
30  format(i10)
    !
    string = obsrv%IDstring
    ! -- Extract node number from string and store it.
    !    If 1st item is not an integer(I4B), it should be a
    !    feature name--deal with it.
    icol = 1
    ! -- get node number
    call urword(string, icol, istart, istop, 1, n, r, iout, inunitobs)
    read (string(istart:istop), 30, iostat=istat) nn
    if (istat == 0) then
      ! -- store uzf node number (NodeNumber)
      obsrv%NodeNumber = nn
    else
      ! Integer can't be read from string; it's presumed to be a boundary
      ! name (already converted to uppercase)
      obsrv%FeatureName = string(istart:istop)
      !obsrv%FeatureName = trim(adjustl(string))
      ! -- Observation may require summing rates from multiple boundaries,
      !    so assign NodeNumber as a value that indicates observation
      !    is for a named boundary or group of boundaries.
      obsrv%NodeNumber = NAMEDBOUNDFLAG
    end if
    !
    ! -- for soil water observation, store depth
    if (obsrv%ObsTypeId == 'WATER-CONTENT') then
      call urword(string, icol, istart, istop, 3, n, r, iout, inunitobs)
      obsdepth = r
      ! -- store observations depth
      obsrv%Obsdepth = obsdepth
    end if
  end subroutine uzf_process_obsID

end submodule UzfModuleObs
