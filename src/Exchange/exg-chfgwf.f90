!> @brief This module contains the ChfGwfExchangeModule Module
!!
!! This module contains a lightweight CHF-GWF exchange class which is
!! primarily based on the underlying and generic SWF-GWF code for connecting
!< a SWF model (either CHF or OLF) with a GWF model.
module ChfGwfExchangeModule

  use KindModule, only: DP, I4B, LGP
  use SimVariablesModule, only: errmsg, iout, model_loc_idx
  use SimModule, only: count_errors, store_error, store_error_filename
  use ListsModule, only: basemodellist, baseexchangelist
  use MemoryHelperModule, only: create_mem_path
  use BaseExchangeModule, only: BaseExchangeType, AddBaseExchangeToList
  use BaseModelModule, only: BaseModelType, GetBaseModelFromList
  use SwfGwfExchangeModule, only: SwfGwfExchangeType
  use GwfModule, only: GwfModelType
  use SwfModule, only: SwfModelType
  use ObsModule, only: obs_cr

  implicit none
  private
  public :: chfgwf_cr

  type, extends(SwfGwfExchangeType) :: ChfGwfExchangeType
  contains
    procedure :: exg_df => chf_gwf_df
    procedure :: source_options
    procedure :: source_dimensions
    procedure :: source_data
  end type ChfGwfExchangeType

contains

  !> @ brief Create CHF GWF exchange
  !!
  !! High level wrapper over the SWF-GWF exchange, which is generic in that
  !! it can be used to couple either a CHF or OLF surface water model to a
  !! GWF model.
  !<
  subroutine chfgwf_cr(filename, name, id, m1_id, m2_id, input_mempath)
    ! modules
    ! dummy
    character(len=*), intent(in) :: filename !< filename for reading
    character(len=*) :: name !< exchange name
    integer(I4B), intent(in) :: id !< id for the exchange
    integer(I4B), intent(in) :: m1_id !< id for model 1
    integer(I4B), intent(in) :: m2_id !< id for model 2
    character(len=*), intent(in) :: input_mempath
    ! local
    type(ChfGwfExchangeType), pointer :: exchange
    class(BaseExchangeType), pointer :: baseexchange

    ! Create a new exchange and add it to the baseexchangelist container
    allocate (exchange)
    baseexchange => exchange
    call AddBaseExchangeToList(baseexchangelist, baseexchange)

    call exchange%initialize(filename, name, 'CHF', id, m1_id, m2_id, &
                             input_mempath)

  end subroutine chfgwf_cr

  !> @ brief Define CHF GWF exchange
  !<
  subroutine chf_gwf_df(this)
    !  modules
    !  dummy
    class(ChfGwfExchangeType) :: this !<  SwfGwfExchangeType
    !  local

    !  log the exchange
    write (iout, '(/a,a)') ' Creating exchange: ', this%name
    !
    !  Ensure models are in same solution
    if (associated(this%swfmodel) .and. associated(this%gwfmodel)) then
      if (this%swfmodel%idsoln /= this%gwfmodel%idsoln) then
        call store_error('Two models are connected in a SWF-GWF '// &
                         'exchange but they are in different solutions. '// &
                         'Models must be in same solution: '// &
                         trim(this%swfmodel%name)//' '// &
                         trim(this%gwfmodel%name))
        call store_error_filename(this%filename)
      end if
    end if

    !  source options
    call this%source_options(iout)

    !  source dimensions
    call this%source_dimensions(iout)

    !  allocate arrays
    call this%allocate_arrays()

    !  source exchange data
    call this%source_data(iout)

    !  Store obs
    ! call this%swf_gwf_df_obs()
    ! if (associated(this%swfmodel1)) then
    !   call this%obs%obs_df(iout, this%name, 'SWF-GWF', this%swfmodel1%dis)
    ! end if

    ! !  validate
    ! call this%validate_exchange()

  end subroutine chf_gwf_df

  !> @ brief Source options
  !<
  subroutine source_options(this, iout)
    !  modules
    use ConstantsModule, only: LENVARNAME, DEM6
    use InputOutputModule, only: getunit, openfile
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    use ExgChfgwfInputModule, only: ExgChfgwfParamFoundType
    use SourceCommonModule, only: filein_fname
    !  dummy
    class(ChfGwfExchangeType) :: this !<  GwfExchangeType
    integer(I4B), intent(in) :: iout
    !  local
    type(ExgChfgwfParamFoundType) :: found

    !  update defaults with idm sourced values
    call mem_set_value(this%ipr_input, 'IPR_INPUT', &
                       this%input_mempath, found%ipr_input)
    call mem_set_value(this%ipr_flow, 'IPR_FLOW', &
                       this%input_mempath, found%ipr_flow)
    call mem_set_value(this%ifixedcond, 'IFIXEDCOND', &
                       this%input_mempath, found%ifixedcond)

    write (iout, '(1x,a)') 'Processing CHF-GWF exchange options'

    if (found%ipr_input) then
      write (iout, '(4x,a)') &
        'The list of exchanges will be printed.'
    end if

    if (found%ipr_flow) then
      write (iout, '(4x,a)') &
        'Exchange flows will be printed to list files.'
    end if

    if (found%ifixedcond) then
      write (iout, '(4x,a)') &
        'Conductance is fixed as product of BEDLEAK and CFACT.'
    end if

    !  enforce 0 or 1 OBS6_FILENAME entries in option block
    ! if (.not. this%is_datacopy) then
    !   if (filein_fname(this%obs%inputFilename, 'OBS6_FILENAME', &
    !                    this%input_mempath, this%filename)) then
    !     this%obs%active = .true.
    !     this%obs%inUnitObs = GetUnit()
    !     call openfile(this%obs%inUnitObs, iout, this%obs%inputFilename, 'OBS')
    !   end if
    ! end if

    write (iout, '(1x,a)') 'End of CHF-GWF exchange options'

  end subroutine source_options

  !> @brief Source dimension from input context
  !<
  subroutine source_dimensions(this, iout)
    !  modules
    use MemoryManagerExtModule, only: mem_set_value
    use ExgChfgwfInputModule, only: ExgChfgwfParamFoundType
    !  dummy
    class(ChfGwfExchangeType) :: this !< instance of exchange object
    integer(I4B), intent(in) :: iout !< for logging
    !  local
    type(ExgChfgwfParamFoundType) :: found

    !  update defaults with idm sourced values
    call mem_set_value(this%nexg, 'NEXG', this%input_mempath, found%nexg)

    write (iout, '(1x,a)') 'PROCESSING EXCHANGE DIMENSIONS'

    if (found%nexg) then
      write (iout, '(4x,a,i0)') 'NEXG = ', this%nexg
    end if

    write (iout, '(1x,a)') 'END OF EXCHANGE DIMENSIONS'

  end subroutine source_dimensions

  !> @brief Source exchange data from input context
  !<
  subroutine source_data(this, iout)
    !  modules
    use MemoryManagerModule, only: mem_setptr
    !  dummy
    class(ChfGwfExchangeType) :: this !< instance of exchange object
    integer(I4B), intent(in) :: iout !< the output file unit
    !  local
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm1
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm2
    real(DP), dimension(:), contiguous, pointer :: bedleak
    real(DP), dimension(:), contiguous, pointer :: cfact
    character(len=20) :: cellstr1, cellstr2
    integer(I4B) :: nerr
    integer(I4B) :: iexg, nodeswf, nodegwf
    !  format
    character(len=*), parameter :: fmtexglabel = "(1x, 3a10, 50(a16))"
    character(len=*), parameter :: fmtexgdata = &
                                   "(5x, a, 1x, a ,50(1pg16.6))"

    call mem_setptr(cellidm1, 'CELLIDM1', this%input_mempath)
    call mem_setptr(cellidm2, 'CELLIDM2', this%input_mempath)
    call mem_setptr(bedleak, 'BEDLEAK', this%input_mempath)
    call mem_setptr(cfact, 'CFACT', this%input_mempath)

    write (iout, '(1x,a)') 'PROCESSING EXCHANGEDATA'

    if (this%ipr_input /= 0) then
      write (iout, fmtexglabel) 'NODEM1', 'NODEM2', 'BEDLEAK', 'CFACT'
    end if

    do iexg = 1, this%nexg

      if (associated(this%model1)) then
        !  Determine user node number
        nodeswf = this%noder(this%model1, cellidm1(:, iexg), iout)
        this%nodeswf(iexg) = nodeswf
      else
        this%nodeswf(iexg) = -1
      end if

      if (associated(this%model2)) then
        !  Determine user node number
        nodegwf = this%noder(this%model2, cellidm2(:, iexg), iout)
        this%nodegwf(iexg) = nodegwf
      else
        this%nodegwf(iexg) = -1
      end if

      !  Read rest of input line
      this%bedleak(iexg) = bedleak(iexg)
      this%cfact(iexg) = cfact(iexg)

      !  Write the data to listing file if requested
      if (this%ipr_input /= 0) then
        cellstr1 = this%cellstr(this%model1, cellidm1(:, iexg), iout)
        cellstr2 = this%cellstr(this%model2, cellidm2(:, iexg), iout)
        write (iout, fmtexgdata) trim(cellstr1), trim(cellstr2), &
          this%bedleak(iexg), this%cfact(iexg)
      end if

      !  Check to see if nodeswf is outside of active domain
      if (associated(this%model1)) then
        if (nodeswf <= 0) then
          cellstr1 = this%cellstr(this%model1, cellidm1(:, iexg), iout)
          write (errmsg, *) &
            trim(adjustl(this%model1%name))// &
            ' Cell is outside active grid domain ('// &
            trim(adjustl(cellstr1))//').'
          call store_error(errmsg)
        end if
      end if

      !  Check to see if nodegwf is outside of active domain
      if (associated(this%model2)) then
        if (nodegwf <= 0) then
          cellstr2 = this%cellstr(this%model2, cellidm2(:, iexg), iout)
          write (errmsg, *) &
            trim(adjustl(this%model2%name))// &
            ' Cell is outside active grid domain ('// &
            trim(adjustl(cellstr2))//').'
          call store_error(errmsg)
        end if
      end if
    end do

    write (iout, '(1x,a)') 'END OF EXCHANGEDATA'

    !  Stop if errors
    nerr = count_errors()
    if (nerr > 0) then
      call store_error('Errors encountered in exchange input file.')
      call store_error_filename(this%filename)
    end if

  end subroutine source_data

end module ChfGwfExchangeModule
