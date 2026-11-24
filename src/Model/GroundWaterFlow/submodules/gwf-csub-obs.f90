!> @brief CSUB package observations
!<
submodule(GwfCsubModule) GwfCsubModuleObs
contains

  ! -- Procedures related to observations (type-bound)
  !> @brief Determine if observations are supported.
  !!
  !! Function to determine if observations are supported by the CSUB package.
  !! Observations are supported by the CSUB package.
  !!
  !<
  module procedure csub_obs_supported
  !
  ! -- initialize variables
  csub_obs_supported = .true.
  end procedure csub_obs_supported

  !> @brief Define the observation types available in the package
  !!
  !! Method to define the observation types available in the CSUB package.
  !!
  !<
  module procedure csub_df_obs
  ! -- local variables
  integer(I4B) :: indx
  !
  ! -- Store obs type and assign procedure pointer
  !    for csub observation type.
  call this%obs%StoreObsType('csub', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for inelastic-csub observation type.
  call this%obs%StoreObsType('inelastic-csub', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for elastic-csub observation type.
  call this%obs%StoreObsType('elastic-csub', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for coarse-csub observation type.
  call this%obs%StoreObsType('coarse-csub', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for csub-cell observation type.
  call this%obs%StoreObsType('csub-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for watercomp-csub observation type.
  call this%obs%StoreObsType('wcomp-csub-cell', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for interbed ske observation type.
  call this%obs%StoreObsType('ske', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for interbed sk observation type.
  call this%obs%StoreObsType('sk', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for ske-cell observation type.
  call this%obs%StoreObsType('ske-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for sk-cell observation type.
  call this%obs%StoreObsType('sk-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for geostatic-stress-cell observation type.
  call this%obs%StoreObsType('gstress-cell', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for effective-stress-cell observation type.
  call this%obs%StoreObsType('estress-cell', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for total-compaction observation type.
  call this%obs%StoreObsType('interbed-compaction', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for inelastic-compaction observation type.
  call this%obs%StoreObsType('inelastic-compaction', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for inelastic-compaction observation type.
  call this%obs%StoreObsType('elastic-compaction', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for coarse-compaction observation type.
  call this%obs%StoreObsType('coarse-compaction', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for inelastic-compaction-cell observation type.
  call this%obs%StoreObsType('inelastic-compaction-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for elastic-compaction-cell observation type.
  call this%obs%StoreObsType('elastic-compaction-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for compaction-cell observation type.
  call this%obs%StoreObsType('compaction-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for interbed thickness observation type.
  call this%obs%StoreObsType('thickness', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for coarse-thickness observation type.
  call this%obs%StoreObsType('coarse-thickness', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for thickness-cell observation type.
  call this%obs%StoreObsType('thickness-cell', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for interbed theta observation type.
  call this%obs%StoreObsType('theta', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for coarse-theta observation type.
  call this%obs%StoreObsType('coarse-theta', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for theta-cell observation type.
  call this%obs%StoreObsType('theta-cell', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for preconstress-cell observation type.
  call this%obs%StoreObsType('preconstress-cell', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for interbed-compaction-pct observation type.
  call this%obs%StoreObsType('interbed-compaction-pct', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-preconstress observation type.
  call this%obs%StoreObsType('delay-preconstress', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-head observation type.
  call this%obs%StoreObsType('delay-head', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-gstress observation type.
  call this%obs%StoreObsType('delay-gstress', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-estress observation type.
  call this%obs%StoreObsType('delay-estress', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-compaction observation type.
  call this%obs%StoreObsType('delay-compaction', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-thickness observation type.
  call this%obs%StoreObsType('delay-thickness', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-theta observation type.
  call this%obs%StoreObsType('delay-theta', .false., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-flowtop observation type.
  call this%obs%StoreObsType('delay-flowtop', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  !
  ! -- Store obs type and assign procedure pointer
  !    for delay-flowbot observation type.
  call this%obs%StoreObsType('delay-flowbot', .true., indx)
  this%obs%obsData(indx)%ProcessIdPtr => csub_process_obsID
  end procedure csub_df_obs

  !> @brief Set the observations for this time step
  !!
  !! Method to set the CSUB package observations for this time step.
  !!
  !<
  module procedure csub_bd_obs
  ! -- local variables
  type(ObserveType), pointer :: obsrv => null()
  integer(I4B) :: i
  integer(I4B) :: j
  integer(I4B) :: n
  integer(I4B) :: idelay
  integer(I4B) :: ncol
  integer(I4B) :: node
  real(DP) :: v
  real(DP) :: r
  real(DP) :: f
  real(DP) :: b0
  !
  ! -- Fill simulated values for all csub observations
  if (this%obs%npakobs > 0) then
    call this%obs%obs_bd_clear()
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      if (obsrv%BndFound) then
        if (obsrv%ObsTypeId == 'SKE' .or. &
            obsrv%ObsTypeId == 'SK' .or. &
            obsrv%ObsTypeId == 'SKE-CELL' .or. &
            obsrv%ObsTypeId == 'SK-CELL' .or. &
            obsrv%ObsTypeId == 'DELAY-HEAD' .or. &
            obsrv%ObsTypeId == 'DELAY-PRECONSTRESS' .or. &
            obsrv%ObsTypeId == 'DELAY-GSTRESS' .or. &
            obsrv%ObsTypeId == 'DELAY-ESTRESS' .or. &
            obsrv%ObsTypeId == 'PRECONSTRESS-CELL') then
          if (this%gwfiss /= 0) then
            call this%obs%SaveOneSimval(obsrv, DNODATA)
          else
            v = DZERO
            do j = 1, obsrv%indxbnds_count
              n = obsrv%indxbnds(j)
              select case (obsrv%ObsTypeId)
              case ('SKE')
                v = this%ske(n)
              case ('SK')
                v = this%sk(n)
              case ('SKE-CELL')
                !
                ! -- add the coarse component
                if (j == 1) then
                  v = this%cg_ske(n)
                else
                  v = this%ske(n)
                end if
              case ('SK-CELL')
                !
                ! -- add the coarse component
                if (j == 1) then
                  v = this%cg_sk(n)
                else
                  v = this%sk(n)
                end if
              case ('DELAY-HEAD', 'DELAY-PRECONSTRESS', &
                    'DELAY-GSTRESS', 'DELAY-ESTRESS')
                if (n > this%ndelaycells) then
                  r = real(n - 1, DP) / real(this%ndelaycells, DP)
                  idelay = int(floor(r)) + 1
                  ncol = n - int(floor(r)) * this%ndelaycells
                else
                  idelay = 1
                  ncol = n
                end if
                select case (obsrv%ObsTypeId)
                case ('DELAY-HEAD')
                  v = this%dbh(ncol, idelay)
                case ('DELAY-PRECONSTRESS')
                  v = this%dbpcs(ncol, idelay)
                case ('DELAY-GSTRESS')
                  v = this%dbgeo(ncol, idelay)
                case ('DELAY-ESTRESS')
                  v = this%dbes(ncol, idelay)
                end select
              case ('PRECONSTRESS-CELL')
                v = this%pcs(n)
              case default
                errmsg = "Unrecognized observation type '"// &
                         trim(obsrv%ObsTypeId)//"'."
                call store_error(errmsg)
              end select
              call this%obs%SaveOneSimval(obsrv, v)
            end do
          end if
        else
          v = DZERO
          do j = 1, obsrv%indxbnds_count
            n = obsrv%indxbnds(j)
            select case (obsrv%ObsTypeId)
            case ('CSUB')
              v = this%storagee(n) + this%storagei(n)
            case ('INELASTIC-CSUB')
              v = this%storagei(n)
            case ('ELASTIC-CSUB')
              v = this%storagee(n)
            case ('COARSE-CSUB')
              v = this%cg_stor(n)
            case ('WCOMP-CSUB-CELL')
              v = this%cell_wcstor(n)
            case ('CSUB-CELL')
              !
              ! -- add the coarse component
              if (j == 1) then
                v = this%cg_stor(n)
              else
                v = this%storagee(n) + this%storagei(n)
              end if
            case ('THETA')
              v = this%theta(n)
            case ('COARSE-THETA')
              v = this%cg_theta(n)
            case ('THETA-CELL')
              !
              ! -- add the coarse component
              if (j == 1) then
                f = this%cg_thick(n) / this%cell_thick(n)
                v = f * this%cg_theta(n)
              else
                node = this%nodelist(n)
                f = this%csub_calc_interbed_thickness(n) / this%cell_thick(node)
                v = f * this%theta(n)
              end if
            case ('GSTRESS-CELL')
              v = this%cg_gs(n)
            case ('ESTRESS-CELL')
              v = this%cg_es(n)
            case ('INTERBED-COMPACTION')
              v = this%tcomp(n)
            case ('INTERBED-COMPACTION-PCT')
              b0 = this%thickini(n)
              if (this%idelay(n) /= 0) then
                b0 = b0 * this%rnb(n)
              end if
              v = DHUNDRED * this%tcomp(n) / b0
            case ('INELASTIC-COMPACTION')
              v = this%tcompi(n)
            case ('ELASTIC-COMPACTION')
              v = this%tcompe(n)
            case ('COARSE-COMPACTION')
              v = this%cg_tcomp(n)
            case ('INELASTIC-COMPACTION-CELL')
              !
              ! -- no coarse inelastic component
              if (j > 1) then
                v = this%tcompi(n)
              end if
            case ('ELASTIC-COMPACTION-CELL')
              !
              ! -- add the coarse component
              if (j == 1) then
                v = this%cg_tcomp(n)
              else
                v = this%tcompe(n)
              end if
            case ('COMPACTION-CELL')
              !
              ! -- add the coarse component
              if (j == 1) then
                v = this%cg_tcomp(n)
              else
                v = this%tcomp(n)
              end if
            case ('THICKNESS')
              idelay = this%idelay(n)
              v = this%thick(n)
              if (idelay /= 0) then
                v = v * this%rnb(n)
              end if
            case ('COARSE-THICKNESS')
              v = this%cg_thick(n)
            case ('THICKNESS-CELL')
              v = this%cell_thick(n)
            case ('DELAY-COMPACTION', 'DELAY-THICKNESS', &
                  'DELAY-THETA')
              if (n > this%ndelaycells) then
                r = real(n, DP) / real(this%ndelaycells, DP)
                idelay = int(floor(r)) + 1
                ncol = mod(n, this%ndelaycells)
              else
                idelay = 1
                ncol = n
              end if
              select case (obsrv%ObsTypeId)
              case ('DELAY-COMPACTION')
                v = this%dbtcomp(ncol, idelay)
              case ('DELAY-THICKNESS')
                v = this%dbdz(ncol, idelay)
              case ('DELAY-THETA')
                v = this%dbtheta(ncol, idelay)
              end select
            case ('DELAY-FLOWTOP')
              idelay = this%idelay(n)
              v = this%dbflowtop(idelay)
            case ('DELAY-FLOWBOT')
              idelay = this%idelay(n)
              v = this%dbflowbot(idelay)
            case default
              errmsg = "Unrecognized observation type: '"// &
                       trim(obsrv%ObsTypeId)//"'."
              call store_error(errmsg)
            end select
            call this%obs%SaveOneSimval(obsrv, v)
          end do
        end if
      else
        call this%obs%SaveOneSimval(obsrv, DNODATA)
      end if
    end do
    !
    ! -- write summary of package error messages
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end if
  end procedure csub_bd_obs

  !> @brief Read and prepare the observations
  !!
  !! Method to read and prepare the observations for the CSUB package.
  !!
  !<
  module procedure csub_rp_obs
  ! -- modules
  use TdisModule, only: kper
  ! -- local variables
  class(ObserveType), pointer :: obsrv => null()
  character(len=LENBOUNDNAME) :: bname
  integer(I4B) :: i
  integer(I4B) :: j
  integer(I4B) :: n
  integer(I4B) :: n2
  integer(I4B) :: idelay
  !
  !  -- return if observations are not supported
  if (.not. this%csub_obs_supported()) then
    return
  end if
  !
  ! -- process each package observation
  !    only done the first stress period since boundaries are fixed
  !    for the simulation
  if (kper == 1) then
    do i = 1, this%obs%npakobs
      obsrv => this%obs%pakobs(i)%obsrv
      !
      ! -- initialize BndFound to .false.
      obsrv%BndFound = .false.
      !
      bname = obsrv%FeatureName
      if (bname /= '') then
        !
        ! -- Observation location(s) is(are) based on a boundary name.
        !    Iterate through all boundaries to identify and store
        !    corresponding index(indices) in bound array.
        if (this%inamedbound /= 0) then
          do j = 1, this%ninterbeds
            if (this%boundname(j) == bname) then
              obsrv%BndFound = .true.
              obsrv%CurrentTimeStepEndValue = DZERO
              call obsrv%AddObsIndex(j)
            end if
          end do
        end if
        if (.not. obsrv%BndFound) then
          write (errmsg, '(7a)') &
            'Boundary "', trim(bname), '" for observation "', &
            trim(obsrv%Name), '" is invalid in package "', &
            trim(this%packName), '"'
          call store_error(errmsg)
        end if
        !
        ! -- one value per cell
      else if (obsrv%ObsTypeId == 'GSTRESS-CELL' .or. &
               obsrv%ObsTypeId == 'ESTRESS-CELL' .or. &
               obsrv%ObsTypeId == 'THICKNESS-CELL' .or. &
               obsrv%ObsTypeId == 'COARSE-CSUB' .or. &
               obsrv%ObsTypeId == 'WCOMP-CSUB-CELL' .or. &
               obsrv%ObsTypeId == 'COARSE-COMPACTION' .or. &
               obsrv%ObsTypeId == 'COARSE-THETA' .or. &
               obsrv%ObsTypeId == 'COARSE-THICKNESS') then
        obsrv%BndFound = .true.
        obsrv%CurrentTimeStepEndValue = DZERO
        call obsrv%AddObsIndex(obsrv%NodeNumber)
      else if (obsrv%ObsTypeId == 'DELAY-PRECONSTRESS' .or. &
               obsrv%ObsTypeId == 'DELAY-HEAD' .or. &
               obsrv%ObsTypeId == 'DELAY-GSTRESS' .or. &
               obsrv%ObsTypeId == 'DELAY-ESTRESS' .or. &
               obsrv%ObsTypeId == 'DELAY-COMPACTION' .or. &
               obsrv%ObsTypeId == 'DELAY-THICKNESS' .or. &
               obsrv%ObsTypeId == 'DELAY-THETA') then
        if (this%ninterbeds > 0) then
          n = obsrv%NodeNumber
          idelay = this%idelay(n)
          if (idelay /= 0) then
            j = (idelay - 1) * this%ndelaycells + 1
            n2 = obsrv%NodeNumber2
            if (n2 < 1 .or. n2 > this%ndelaycells) then
              write (errmsg, '(a,2(1x,a),1x,i0,1x,a,i0,a)') &
                trim(adjustl(obsrv%ObsTypeId)), 'interbed cell must be ', &
                'greater than 0 and less than or equal to', this%ndelaycells, &
                '(specified value is ', n2, ').'
              call store_error(errmsg)
            else
              j = (idelay - 1) * this%ndelaycells + n2
            end if
            obsrv%BndFound = .true.
            call obsrv%AddObsIndex(j)
          end if
        end if
        !
        ! -- interbed value
      else if (obsrv%ObsTypeId == 'CSUB' .or. &
               obsrv%ObsTypeId == 'INELASTIC-CSUB' .or. &
               obsrv%ObsTypeId == 'ELASTIC-CSUB' .or. &
               obsrv%ObsTypeId == 'SK' .or. &
               obsrv%ObsTypeId == 'SKE' .or. &
               obsrv%ObsTypeId == 'THICKNESS' .or. &
               obsrv%ObsTypeId == 'THETA' .or. &
               obsrv%ObsTypeId == 'INTERBED-COMPACTION' .or. &
               obsrv%ObsTypeId == 'INELASTIC-COMPACTION' .or. &
               obsrv%ObsTypeId == 'ELASTIC-COMPACTION' .or. &
               obsrv%ObsTypeId == 'INTERBED-COMPACTION-PCT') then
        if (this%ninterbeds > 0) then
          j = obsrv%NodeNumber
          if (j < 1 .or. j > this%ninterbeds) then
            write (errmsg, '(a,2(1x,a),1x,i0,1x,a,i0,a)') &
              trim(adjustl(obsrv%ObsTypeId)), 'interbed cell must be greater', &
              'than 0 and less than or equal to', this%ninterbeds, &
              '(specified value is ', j, ').'
            call store_error(errmsg)
          else
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(j)
          end if
        end if
      else if (obsrv%ObsTypeId == 'DELAY-FLOWTOP' .or. &
               obsrv%ObsTypeId == 'DELAY-FLOWBOT') then
        if (this%ninterbeds > 0) then
          j = obsrv%NodeNumber
          if (j < 1 .or. j > this%ninterbeds) then
            write (errmsg, '(a,2(1x,a),1x,i0,1x,a,i0,a)') &
              trim(adjustl(obsrv%ObsTypeId)), &
              'interbed cell must be greater ', &
              'than 0 and less than or equal to', this%ninterbeds, &
              '(specified value is ', j, ').'
            call store_error(errmsg)
          end if
          idelay = this%idelay(j)
          if (idelay /= 0) then
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(j)
          end if
        end if
      else
        !
        ! -- Accumulate values in a single cell
        ! -- Observation location is a single node number
        ! -- save node number in first position
        if (obsrv%ObsTypeId == 'CSUB-CELL' .or. &
            obsrv%ObsTypeId == 'SKE-CELL' .or. &
            obsrv%ObsTypeId == 'SK-CELL' .or. &
            obsrv%ObsTypeId == 'THETA-CELL' .or. &
            obsrv%ObsTypeId == 'INELASTIC-COMPACTION-CELL' .or. &
            obsrv%ObsTypeId == 'ELASTIC-COMPACTION-CELL' .or. &
            obsrv%ObsTypeId == 'COMPACTION-CELL') then
          if (.NOT. obsrv%BndFound) then
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(obsrv%NodeNumber)
          end if
        end if
        jloop: do j = 1, this%ninterbeds
          if (this%nodelist(j) == obsrv%NodeNumber) then
            obsrv%BndFound = .true.
            obsrv%CurrentTimeStepEndValue = DZERO
            call obsrv%AddObsIndex(j)
          end if
        end do jloop
      end if
    end do
    !
    ! -- evaluate if there are any observation errors
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end if
  end procedure csub_rp_obs

  ! -- Procedures related to observations (NOT type-bound)
  !> @brief Process the observation IDs for the package
  !!
  !! Method to process the observation IDs for the CSUB package. This
  !! procedure is pointed to by ObsDataType%ProcesssIdPtr. It processes the
  !! ID string of an observation definition for csub-package observations.
  !!
  !<
  module procedure csub_process_obsID
  ! -- local variables
  integer(I4B) :: nn1
  integer(I4B) :: nn2
  integer(I4B) :: icol, istart, istop
  character(len=LINELENGTH) :: string
  character(len=LENBOUNDNAME) :: bndname
  logical(LGP) :: flag_string
  logical(LGP) :: flag_idcellno
  logical(LGP) :: flag_error
  !
  ! -- initialize variables
  string = obsrv%IDstring
  flag_string = .TRUE.
  flag_idcellno = .FALSE.
  flag_error = .FALSE.
  if (obsrv%ObsTypeId(1:5) == "DELAY" .AND. &
      obsrv%ObsTypeId(1:10) /= "DELAY-FLOW") then
    flag_idcellno = .TRUE.
  end if
  !
  ! -- Extract reach number from string and store it.
  !    If 1st item is not an integer(I4B), it should be a
  !    boundary name--deal with it.
  icol = 1
  !
  ! -- get icsubno number or boundary name
  if (obsrv%ObsTypeId == 'CSUB' .or. &
      obsrv%ObsTypeId == 'INELASTIC-CSUB' .or. &
      obsrv%ObsTypeId == 'ELASTIC-CSUB' .or. &
      obsrv%ObsTypeId == 'SK' .or. &
      obsrv%ObsTypeId == 'SKE' .or. &
      obsrv%ObsTypeId == 'THETA' .or. &
      obsrv%ObsTypeId == 'THICKNESS' .or. &
      obsrv%ObsTypeId == 'INTERBED-COMPACTION' .or. &
      obsrv%ObsTypeId == 'INTERBED-COMPACTION-PCT' .or. &
      obsrv%ObsTypeId == 'INELASTIC-COMPACTION' .or. &
      obsrv%ObsTypeId == 'ELASTIC-COMPACTION' .or. &
      obsrv%ObsTypeId == 'DELAY-HEAD' .or. &
      obsrv%ObsTypeId == 'DELAY-GSTRESS' .or. &
      obsrv%ObsTypeId == 'DELAY-ESTRESS' .or. &
      obsrv%ObsTypeId == 'DELAY-PRECONSTRESS' .or. &
      obsrv%ObsTypeId == 'DELAY-COMPACTION' .or. &
      obsrv%ObsTypeId == 'DELAY-THICKNESS' .or. &
      obsrv%ObsTypeId == 'DELAY-THETA' .or. &
      obsrv%ObsTypeId == 'DELAY-FLOWTOP' .or. &
      obsrv%ObsTypeId == 'DELAY-FLOWBOT') then
    call extract_idnum_or_bndname(string, icol, istart, istop, nn1, bndname)
    ! read cellid
  else
    nn1 = dis%noder_from_string(icol, istart, istop, inunitobs, &
                                iout, string, flag_string)
  end if
  ! boundnames are not allowed for these observation types
  if (obsrv%ObsTypeId == 'SK' .or. &
      obsrv%ObsTypeId == 'SKE' .or. &
      obsrv%ObsTypeId == 'THETA' .or. &
      obsrv%ObsTypeId == 'THICKNESS' .or. &
      obsrv%ObsTypeId == 'INTERBED-COMPACTION' .or. &
      obsrv%ObsTypeId == 'INELASTIC-COMPACTION' .or. &
      obsrv%ObsTypeId == 'ELASTIC-COMPACTION' .or. &
      obsrv%ObsTypeId == 'DELAY-HEAD' .or. &
      obsrv%ObsTypeId == 'DELAY-GSTRESS' .or. &
      obsrv%ObsTypeId == 'DELAY-ESTRESS' .or. &
      obsrv%ObsTypeId == 'DELAY-PRECONSTRESS' .or. &
      obsrv%ObsTypeId == 'DELAY-COMPACTION' .or. &
      obsrv%ObsTypeId == 'DELAY-THICKNESS' .or. &
      obsrv%ObsTypeId == 'DELAY-THETA') then
    if (nn1 == NAMEDBOUNDFLAG) then
      write (errmsg, '(5a)') &
        "BOUNDNAME ('", trim(adjustl(bndname)), &
        "') not allowed for CSUB observation type '", &
        trim(adjustl(obsrv%ObsTypeId)), "'."
      call store_error(errmsg)
      flag_error = .TRUE.
    end if
    ! boundnames are allowed for these observation types
  else if (obsrv%ObsTypeId == 'CSUB' .or. &
           obsrv%ObsTypeId == 'INELASTIC-CSUB' .or. &
           obsrv%ObsTypeId == 'ELASTIC-CSUB' .or. &
           !  obsrv%ObsTypeId == 'INTERBED-COMPACTION' .or. &
           !  obsrv%ObsTypeId == 'INELASTIC-COMPACTION' .or. &
           !  obsrv%ObsTypeId == 'ELASTIC-COMPACTION' .or. &
           obsrv%ObsTypeId == 'DELAY-FLOWTOP' .or. &
           obsrv%ObsTypeId == 'DELAY-FLOWBOT') then
    if (nn1 == NAMEDBOUNDFLAG) then
      obsrv%FeatureName = bndname
    end if
  end if
  ! read idcellno for delay observations
  if (flag_idcellno .EQV. .TRUE. .AND. flag_error .EQV. .FALSE.) then
    if (nn1 /= NAMEDBOUNDFLAG) then
      call extract_idnum_or_bndname(string, icol, istart, istop, nn2, bndname)
      if (nn2 == NAMEDBOUNDFLAG) then
        write (errmsg, '(5a)') &
          "BOUNDNAME ('", trim(adjustl(bndname)), &
          "')  not allowed for CSUB observation type '", &
          trim(adjustl(obsrv%ObsTypeId)), "' idcellno."
        call store_error(errmsg)
      else
        obsrv%NodeNumber2 = nn2
      end if
    end if
  end if
  !
  ! -- store reach number (NodeNumber)
  obsrv%NodeNumber = nn1
  end procedure csub_process_obsID

end submodule GwfCsubModuleObs
