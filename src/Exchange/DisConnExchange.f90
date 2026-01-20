module DisConnExchangeModule
  use KindModule, only: I4B, DP, LGP
  use SimVariablesModule, only: errmsg
  use ConstantsModule, only: LENAUXNAME, LENBOUNDNAME, LINELENGTH
  use SimModule, only: store_error, count_errors, store_error_filename
  use CharacterStringModule
  use ListModule, only: ListType
  use MemoryManagerModule, only: mem_allocate, mem_reallocate
  use NumericalModelModule, only: NumericalModelType
  use VirtualModelModule, only: VirtualModelType
  use NumericalExchangeModule, only: NumericalExchangeType
  implicit none

  private
  public :: DisConnExchangeType
  public :: CastAsDisConnExchangeClass, AddDisConnExchangeToList, &
            GetDisConnExchangeFromList

  !> Exchange based on connection between discretizations of DisBaseType.
  !! The data specifies the connections, similar to the information stored
  !! in the connections object: DisBaseType%con
  !<
  type, extends(NumericalExchangeType) :: DisConnExchangeType
    character(len=LINELENGTH), pointer :: filename => null() !< name of the input file

    class(NumericalModelType), pointer :: model1 => null() !< model 1
    class(NumericalModelType), pointer :: model2 => null() !< model 2
    class(VirtualModelType), pointer :: v_model1 => null() !< virtual model 1
    class(VirtualModelType), pointer :: v_model2 => null() !< virtual model 2
    logical(LGP) :: is_datacopy !< when true, this exchange is just a data copy on another process and
                                !! not responsible for controlling movers, observations, ...

    integer(I4B), pointer :: nexg => null() !< number of exchanges
    integer(I4B), dimension(:), pointer, contiguous :: nodem1 => null() !< node numbers in model 1
    integer(I4B), dimension(:), pointer, contiguous :: nodem2 => null() !< node numbers in model 2
    integer(I4B), dimension(:), pointer, contiguous :: ihc => null() !< horizontal connection indicator array, size: nexg
    real(DP), dimension(:), pointer, contiguous :: cl1 => null() !< connection length 1, size: nexg
    real(DP), dimension(:), pointer, contiguous :: cl2 => null() !< connection length 2, size: nexg
    real(DP), dimension(:), pointer, contiguous :: hwva => null() !< horizontal widths, vertical flow areas, size: nexg
    integer(I4B), pointer :: naux => null() !< number of auxiliary variables
    character(len=LENBOUNDNAME), dimension(:), &
      pointer, contiguous :: boundname => null() !< boundnames

    character(len=LENAUXNAME), dimension(:), &
      pointer, contiguous :: auxname => null() !< vector of auxname
    type(CharacterStringType), dimension(:), pointer, &
      contiguous :: auxname_cst => null() !< copy of vector auxname that can be stored in memory manager
    real(DP), dimension(:, :), pointer, contiguous :: auxvar => null() !< array of auxiliary variable values
    integer(I4B), pointer :: ianglex => null() !< flag indicating anglex was read, if read, ianglex is index in auxvar
    integer(I4B), pointer :: icdist => null() !< flag indicating cdist was read, if read, icdist is index in auxvar
    integer(I4B), pointer :: iprpak => null() !< print input flag
    integer(I4B), pointer :: iprflow => null() !< print flag for cell by cell flows
    integer(I4B), pointer :: ipakcb => null() !< save flag for cell by cell flows
    integer(I4B), pointer :: inamedbound => null() !< flag to read boundnames

    integer(I4B), pointer :: ixt3d => null() !< flag indicating if XT3D should be applied on the interface: 0 = off, 1 = lhs, 2 = rhs
    logical(LGP), pointer :: dev_ifmod_on !< development option, forces interface model for this exchange

  contains

    procedure :: allocate_scalars
    procedure :: allocate_arrays
    procedure :: disconnex_da
    procedure :: use_interface_model

    ! protected
    procedure, pass(this) :: source_options
    procedure, pass(this) :: source_dimensions
    procedure, pass(this) :: source_data
    procedure, pass(this) :: noder
    procedure, pass(this) :: cellstr

  end type DisConnExchangeType

  !> @ brief DisConnExchangeFoundType
  !!
  !!  This type is used to simplify the tracking of common parameters
  !!  that are sourced from the input context.
  !<
  type DisConnExchangeFoundType
    logical :: naux = .false.
    logical :: ipakcb = .false.
    logical :: iprpak = .false.
    logical :: iprflow = .false.
    logical :: boundnames = .false.
    logical :: auxiliary = .false.
    logical :: dev_ifmod_on = .false.
    logical :: nexg = .false.
  end type DisConnExchangeFoundType

contains

  !> @brief Source options from input context
  !<
  subroutine source_options(this, iout)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use ArrayHandlersModule, only: ifind
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    integer(I4B), intent(in) :: iout !< for logging
    ! -- local
    type(DisConnExchangeFoundType) :: found
    integer(I4B) :: ival, n
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%naux, 'NAUX', this%input_mempath, found%naux)
    call mem_set_value(this%ipakcb, 'IPAKCB', this%input_mempath, found%ipakcb)
    call mem_set_value(this%iprpak, 'IPRPAK', this%input_mempath, found%iprpak)
    call mem_set_value(this%iprflow, 'IPRFLOW', this%input_mempath, found%iprflow)
    call mem_set_value(this%inamedbound, 'BOUNDNAMES', this%input_mempath, &
                       found%boundnames)
    call mem_set_value(this%dev_ifmod_on, 'DEV_IFMOD_ON', this%input_mempath, &
                       found%dev_ifmod_on)
    !
    ! -- reallocate aux arrays if aux variables provided
    if (found%naux .and. this%naux > 0) then
      call mem_reallocate(this%auxname, LENAUXNAME, this%naux, &
                          'AUXNAME', this%memoryPath)
      call mem_reallocate(this%auxname_cst, LENAUXNAME, this%naux, &
                          'AUXNAME_CST', this%memoryPath)
      call mem_set_value(this%auxname_cst, 'AUXILIARY', this%input_mempath, &
                         found%auxiliary)
      !
      do n = 1, this%naux
        this%auxname(n) = this%auxname_cst(n)
      end do
      !
      ! -- If ANGLDEGX is an auxiliary variable, then anisotropy can be
      !    used in either model.  Store ANGLDEGX position in this%ianglex
      ival = ifind(this%auxname, 'ANGLDEGX')
      if (ival > 0) then
        this%ianglex = ival
      end if
      !
      ival = ifind(this%auxname, 'CDIST')
      if (ival > 0) then
        this%icdist = ival
      end if
    end if
    !
    if (found%ipakcb) then
      this%ipakcb = -1
      write (iout, '(4x,a)') &
        'EXCHANGE FLOWS WILL BE SAVED TO BINARY BUDGET FILES.'
    end if
    !
    if (found%iprpak) then
      write (iout, '(4x,a)') &
        'THE LIST OF EXCHANGES WILL BE PRINTED.'
    end if
    !
    if (found%iprflow) then
      write (iout, '(4x,a)') &
        'EXCHANGE FLOWS WILL BE PRINTED TO LIST FILES.'
    end if
    !
    if (found%boundnames) then
      write (iout, '(4x,a)') 'EXCHANGE BOUNDARIES HAVE NAMES IN LAST COLUMN'
    end if
    !
    if (found%dev_ifmod_on) then
      write (iout, '(4x,2a)') 'Interface model coupling approach manually &
        &activated for ', trim(this%name)
    end if
  end subroutine source_options

  !> @brief Source dimension from input context
  !<
  subroutine source_dimensions(this, iout)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    integer(I4B), intent(in) :: iout !< for logging
    ! -- local
    type(DisConnExchangeFoundType) :: found
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%nexg, 'NEXG', this%input_mempath, found%nexg)
    !
    write (iout, '(1x,a)') 'PROCESSING EXCHANGE DIMENSIONS'
    !
    if (found%nexg) then
      write (iout, '(4x,a,i0)') 'NEXG = ', this%nexg
    end if
    !
    write (iout, '(1x,a)') 'END OF EXCHANGE DIMENSIONS'
  end subroutine source_dimensions

  !> @brief Returns reduced node number from user
  !< specified cell id.
  function noder(this, model, cellid, iout)
    ! -- modules
    use GeomUtilModule, only: get_node
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    class(NumericalModelType), pointer, intent(in) :: model
    integer(I4B), dimension(:), intent(in) :: cellid
    integer(I4B), intent(in) :: iout !< the output file unit
    integer(I4B) :: noder, node
    !
    if (model%dis%ndim == 1) then
      node = cellid(1)
    elseif (model%dis%ndim == 2) then
      node = get_node(cellid(1), 1, cellid(2), &
                      model%dis%mshape(1), 1, &
                      model%dis%mshape(2))
    else
      node = get_node(cellid(1), cellid(2), cellid(3), &
                      model%dis%mshape(1), &
                      model%dis%mshape(2), &
                      model%dis%mshape(3))
    end if
    noder = model%dis%get_nodenumber(node, 0)
  end function noder

  !> @brief
  !<
  function cellstr(this, ndim, cellid, iout)
    ! -- modules
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    integer(I4B) :: ndim !< model DIS dimension
    integer(I4B), dimension(:), intent(in) :: cellid
    integer(I4B), intent(in) :: iout !< the output file unit
    character(len=20) :: cellstr
    character(len=*), parameter :: fmtndim1 = &
                                   "('(',i0,')')"
    character(len=*), parameter :: fmtndim2 = &
                                   "('(',i0,',',i0,')')"
    character(len=*), parameter :: fmtndim3 = &
                                   "('(',i0,',',i0,',',i0,')')"
    !
    cellstr = ''
    !
    select case (ndim)
    case (1)
      write (cellstr, fmtndim1) cellid(1)
    case (2)
      write (cellstr, fmtndim2) cellid(1), cellid(2)
    case (3)
      write (cellstr, fmtndim3) cellid(1), cellid(2), cellid(3)
    case default
    end select
  end function cellstr

  !> @brief Source exchange data from input context
  !<
  subroutine source_data(this, iout)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    integer(I4B), intent(in) :: iout !< the output file unit
    ! -- local
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm1
    integer(I4B), dimension(:, :), contiguous, pointer :: cellidm2
    integer(I4B), dimension(:), contiguous, pointer :: ihc
    real(DP), dimension(:), contiguous, pointer :: cl1
    real(DP), dimension(:), contiguous, pointer :: cl2
    real(DP), dimension(:), contiguous, pointer :: hwva
    real(DP), dimension(:, :), contiguous, pointer :: auxvar
    type(CharacterStringType), dimension(:), contiguous, pointer :: boundname
    integer(I4B) :: ndim1, ndim2
    character(len=20) :: cellstr1, cellstr2
    character(len=2) :: cnfloat
    integer(I4B) :: nerr, iaux
    integer(I4B) :: iexg, nodem1, nodem2
    ! -- format
    character(len=*), parameter :: fmtexglabel = "(1x, 3a10, 50(a16))"
    character(len=*), parameter :: fmtexgdata = &
                                   "(5x, a, 1x, a ,I10, 50(1pg16.6))"
    character(len=40) :: fmtexgdata2
    !
    call mem_setptr(cellidm1, 'CELLIDM1', this%input_mempath)
    call mem_setptr(cellidm2, 'CELLIDM2', this%input_mempath)
    call mem_setptr(ihc, 'IHC', this%input_mempath)
    call mem_setptr(cl1, 'CL1', this%input_mempath)
    call mem_setptr(cl2, 'CL2', this%input_mempath)
    call mem_setptr(hwva, 'HWVA', this%input_mempath)
    call mem_setptr(auxvar, 'AUXVAR', this%input_mempath)
    call mem_setptr(boundname, 'BOUNDNAME', this%input_mempath)
    ndim1 = size(cellidm1, dim=1)
    ndim2 = size(cellidm2, dim=1)
    !
    write (iout, '(1x,a)') 'PROCESSING EXCHANGEDATA'
    !
    if (this%iprpak /= 0) then
      if (this%inamedbound == 0) then
        write (iout, fmtexglabel) 'NODEM1', 'NODEM2', 'IHC', &
          'CL1', 'CL2', 'HWVA', (adjustr(this%auxname(iaux)), &
                                 iaux=1, this%naux)
      else
        write (iout, fmtexglabel) 'NODEM1', 'NODEM2', 'IHC', 'CL1', 'CL2', &
          'HWVA', (adjustr(this%auxname(iaux)), iaux=1, this%naux), &
          ' BOUNDNAME      '
        ! Define format suitable for writing input data,
        ! any auxiliary variables, and boundname.
        write (cnfloat, '(i0)') 3 + this%naux
        fmtexgdata2 = '(5x, a, 1x, a, i10, '//trim(cnfloat)// &
                      '(1pg16.6), 1x, a)'
      end if
    end if
    !
    do iexg = 1, this%nexg
      !
      if (associated(this%model1)) then
        !
        ! -- Determine reduced node number
        nodem1 = this%noder(this%model1, cellidm1(:, iexg), iout)
        this%nodem1(iexg) = nodem1
        !
      else
        this%nodem1(iexg) = -1
      end if
      !
      if (associated(this%model2)) then
        !
        ! -- Determine reduced node number
        nodem2 = this%noder(this%model2, cellidm2(:, iexg), iout)
        this%nodem2(iexg) = nodem2
        !
      else
        this%nodem2(iexg) = -1
      end if
      !
      ! -- Read rest of input line
      this%ihc(iexg) = ihc(iexg)
      this%cl1(iexg) = cl1(iexg)
      this%cl2(iexg) = cl2(iexg)
      this%hwva(iexg) = hwva(iexg)
      do iaux = 1, this%naux
        this%auxvar(iaux, iexg) = auxvar(iaux, iexg)
      end do
      if (this%inamedbound == 1) then
        this%boundname(iexg) = boundname(iexg)
      end if
      !
      ! -- Write the data to listing file if requested
      if (this%iprpak /= 0) then
        cellstr1 = this%cellstr(ndim1, cellidm1(:, iexg), iout)
        cellstr2 = this%cellstr(ndim2, cellidm2(:, iexg), iout)
        if (this%inamedbound == 0) then
          write (iout, fmtexgdata) trim(cellstr1), trim(cellstr2), &
            this%ihc(iexg), this%cl1(iexg), this%cl2(iexg), &
            this%hwva(iexg), &
            (this%auxvar(iaux, iexg), iaux=1, this%naux)
        else
          write (iout, fmtexgdata2) trim(cellstr1), trim(cellstr2), &
            this%ihc(iexg), this%cl1(iexg), this%cl2(iexg), &
            this%hwva(iexg), &
            (this%auxvar(iaux, iexg), iaux=1, this%naux), &
            trim(this%boundname(iexg))
        end if
      end if
      !
      ! -- Check to see if nodem1 is outside of active domain
      if (associated(this%model1)) then
        if (nodem1 <= 0) then
          cellstr1 = this%cellstr(ndim1, cellidm1(:, iexg), iout)
          write (errmsg, *) &
            trim(adjustl(this%model1%name))// &
            ' Cell is outside active grid domain ('// &
            trim(adjustl(cellstr1))//').'
          call store_error(errmsg)
        end if
      end if
      !
      ! -- Check to see if nodem2 is outside of active domain
      if (associated(this%model2)) then
        if (nodem2 <= 0) then
          cellstr2 = this%cellstr(ndim2, cellidm2(:, iexg), iout)
          write (errmsg, *) &
            trim(adjustl(this%model2%name))// &
            ' Cell is outside active grid domain ('// &
            trim(adjustl(cellstr2))//').'
          call store_error(errmsg)
        end if
      end if
    end do
    !
    write (iout, '(1x,a)') 'END OF EXCHANGEDATA'
    !
    ! -- Stop if errors
    nerr = count_errors()
    if (nerr > 0) then
      call store_error('Errors encountered in exchange input file.')
      call store_error_filename(this%filename)
    end if
  end subroutine source_data

  !> @brief Allocate scalars and initialize to defaults
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    !
    allocate (this%filename)
    this%filename = ''
    !
    call mem_allocate(this%nexg, 'NEXG', this%memoryPath)
    call mem_allocate(this%naux, 'NAUX', this%memoryPath)
    call mem_allocate(this%ianglex, 'IANGLEX', this%memoryPath)
    call mem_allocate(this%icdist, 'ICDIST', this%memoryPath)
    call mem_allocate(this%ixt3d, 'IXT3D', this%memoryPath)
    call mem_allocate(this%iprpak, 'IPRPAK', this%memoryPath)
    call mem_allocate(this%iprflow, 'IPRFLOW', this%memoryPath)
    call mem_allocate(this%ipakcb, 'IPAKCB', this%memoryPath)
    call mem_allocate(this%inamedbound, 'INAMEDBOUND', this%memoryPath)
    call mem_allocate(this%dev_ifmod_on, 'DEV_IFMOD_ON', this%memoryPath)

    call mem_allocate(this%auxname, LENAUXNAME, 0, &
                      'AUXNAME', this%memoryPath)
    call mem_allocate(this%auxname_cst, LENAUXNAME, 0, &
                      'AUXNAME_CST', this%memoryPath)
    !
    this%nexg = 0
    this%naux = 0
    this%ianglex = 0
    this%icdist = 0
    this%ixt3d = 0
    this%iprpak = 0
    this%iprflow = 0
    this%ipakcb = 0
    this%inamedbound = 0
    !
    this%dev_ifmod_on = .false.
  end subroutine allocate_scalars

  !> @brief Allocate array data, using the number of
  !! connected nodes @param nexg
  !<
  subroutine allocate_arrays(this)
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    !
    call mem_allocate(this%nodem1, this%nexg, 'NODEM1', this%memoryPath)
    call mem_allocate(this%nodem2, this%nexg, 'NODEM2', this%memoryPath)
    call mem_allocate(this%ihc, this%nexg, 'IHC', this%memoryPath)
    call mem_allocate(this%cl1, this%nexg, 'CL1', this%memoryPath)
    call mem_allocate(this%cl2, this%nexg, 'CL2', this%memoryPath)
    call mem_allocate(this%hwva, this%nexg, 'HWVA', this%memoryPath)
    ! NB: auxname array is allocated while parsing
    call mem_allocate(this%auxvar, this%naux, this%nexg, &
                      'AUXVAR', this%memoryPath)
    !
    ! allocate boundname
    if (this%inamedbound == 1) then
      allocate (this%boundname(this%nexg))
    else
      allocate (this%boundname(1))
    end if
    this%boundname(:) = ''
  end subroutine allocate_arrays

  !> @brief Should interface model be used to handle these
  !! exchanges, to be overridden for inheriting types
  !<
  function use_interface_model(this) result(use_im)
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    ! -- return
    logical(LGP) :: use_im !< flag whether interface model should be used
                          !! for this exchange instead
    !
    ! use im when one of the models is not local
    use_im = .not. (this%v_model1%is_local .and. this%v_model2%is_local)
  end function use_interface_model

  !> @brief Clean up all scalars and arrays
  !<
  subroutine disconnex_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(DisConnExchangeType) :: this !< instance of exchange object
    !
    ! arrays
    call mem_deallocate(this%nodem1)
    call mem_deallocate(this%nodem2)
    call mem_deallocate(this%ihc)
    call mem_deallocate(this%cl1)
    call mem_deallocate(this%cl2)
    call mem_deallocate(this%hwva)
    call mem_deallocate(this%auxvar)
    !
    deallocate (this%boundname)
    !
    ! scalars
    call mem_deallocate(this%nexg)
    call mem_deallocate(this%naux)
    call mem_deallocate(this%auxname, 'AUXNAME', this%memoryPath)
    call mem_deallocate(this%auxname_cst, 'AUXNAME_CST', this%memoryPath)
    call mem_deallocate(this%ianglex)
    call mem_deallocate(this%icdist)
    call mem_deallocate(this%ixt3d)
    call mem_deallocate(this%iprpak)
    call mem_deallocate(this%iprflow)
    call mem_deallocate(this%ipakcb)
    call mem_deallocate(this%inamedbound)
    call mem_deallocate(this%dev_ifmod_on)
  end subroutine disconnex_da

  function CastAsDisConnExchangeClass(obj) result(res)
    implicit none
    ! -- dummy
    class(*), pointer, intent(inout) :: obj
    ! -- return
    class(DisConnExchangeType), pointer :: res
    !
    res => null()
    if (.not. associated(obj)) return
    !
    select type (obj)
    class is (DisConnExchangeType)
      res => obj
    end select
  end function CastAsDisConnExchangeClass

  subroutine AddDisConnExchangeToList(list, exchange)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    class(DisConnExchangeType), pointer, intent(in) :: exchange
    ! -- local
    class(*), pointer :: obj
    !
    obj => exchange
    call list%Add(obj)
  end subroutine AddDisConnExchangeToList

  function GetDisConnExchangeFromList(list, idx) result(res)
    implicit none
    ! -- dummy
    type(ListType), intent(inout) :: list
    integer(I4B), intent(in) :: idx
    ! -- return
    class(DisConnExchangeType), pointer :: res
    ! -- local
    class(*), pointer :: obj
    !
    obj => list%GetItem(idx)
    res => CastAsDisConnExchangeClass(obj)
  end function GetDisConnExchangeFromList

end module DisConnExchangeModule
