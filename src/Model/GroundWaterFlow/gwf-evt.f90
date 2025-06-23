module EvtModule
  !
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DONE, LENFTYPE, LENPACKAGENAME, MAXCHARLEN, &
                             IWETLAKE
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use SimModule, only: store_error, store_error_filename, count_errors
  use SimVariablesModule, only: errmsg
  use ObsModule, only: DefaultObsIdProcessor
  use CharacterStringModule, only: CharacterStringType
  use MatrixBaseModule
  use GeomUtilModule, only: get_node
  !
  implicit none
  !
  private
  public :: evt_create
  !
  character(len=LENFTYPE) :: ftype = 'EVT'
  character(len=LENPACKAGENAME) :: text = '             EVT'
  character(len=LENPACKAGENAME) :: texta = '           EVTA'
  !
  type, extends(BndExtType) :: EvtType
    ! -- logicals
    logical, pointer, private :: segsdefined
    logical, pointer, private :: fixed_cell
    logical, pointer, private :: read_as_arrays
    logical, pointer, private :: surfratespecified
    ! -- integers
    integer(I4B), pointer, private :: nseg => null() !< number of ET segments
    ! -- arrays
    real(DP), dimension(:), pointer, contiguous :: surface => null() !< elevation of the ET surface
    real(DP), dimension(:), pointer, contiguous :: rate => null() !< maximum ET flux rate
    real(DP), dimension(:), pointer, contiguous :: depth => null() !< ET extinction depth
    real(DP), dimension(:, :), pointer, contiguous :: pxdp => null() !< proportion of ET extinction depth at bottom of segment
    real(DP), dimension(:, :), pointer, contiguous :: petm => null() !< proportion of max ET flux rate at bottom of segment
    real(DP), dimension(:), pointer, contiguous :: petm0 => null() !< proportion of max ET flux rate that will apply when head is at or above ET surface
    integer(I4B), dimension(:), pointer, contiguous :: nodesontop => null()
  contains
    procedure :: evt_allocate_scalars
    procedure :: allocate_arrays => evt_allocate_arrays
    procedure :: source_options => evt_source_options
    procedure :: source_dimensions => evt_source_dimensions
    procedure :: evt_log_options
    procedure :: read_initial_attr => evt_read_initial_attr
    procedure :: bnd_rp => evt_rp
    procedure :: set_nodesontop
    procedure :: bnd_cf => evt_cf
    procedure :: bnd_fc => evt_fc
    procedure :: bnd_da => evt_da
    procedure :: define_listlabel => evt_define_listlabel
    procedure :: bound_value => evt_bound_value
    procedure, private :: default_nodelist
    procedure, private :: check_pxdp
    ! -- for observations
    procedure, public :: bnd_obs_supported => evt_obs_supported
    procedure, public :: bnd_df_obs => evt_df_obs
  end type EvtType

contains

  !> @brief Create a new Evapotranspiration Segments Package and point pakobj
  !! to the new package
  !<
  subroutine evt_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath)
    ! -- dummy
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    ! -- local
    type(EvtType), pointer :: evtobj
    !
    ! -- allocate evt object and scalar variables
    allocate (evtobj)
    packobj => evtobj
    !
    ! -- create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    !
    ! -- allocate scalars
    call evtobj%evt_allocate_scalars()
    !
    ! -- initialize package
    call packobj%pack_initialize()

    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF')
  end subroutine evt_create

  !> @brief Allocate package scalar members
  !<
  subroutine evt_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(EvtType), intent(inout) :: this
    !
    ! -- call standard BndType allocate scalars
    call this%BndExtType%allocate_scalars()
    !
    ! -- allocate the object and assign values to object variables
    call mem_allocate(this%nseg, 'NSEG', this%memoryPath)
    !
    ! -- allocate internal members
    allocate (this%segsdefined)
    allocate (this%fixed_cell)
    allocate (this%read_as_arrays)
    allocate (this%surfratespecified)
    !
    ! -- Set values
    this%nseg = 1
    this%segsdefined = .true.
    this%fixed_cell = .false.
    this%read_as_arrays = .false.
    this%surfratespecified = .false.
  end subroutine evt_allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine evt_allocate_arrays(this, nodelist, auxvar)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr, mem_checkin
    ! -- dummy
    class(EvtType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    !
    ! -- call standard BndType allocate scalars
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    !
    ! -- set EVT input context pointers
    call mem_setptr(this%surface, 'SURFACE', this%input_mempath)
    call mem_setptr(this%rate, 'RATE', this%input_mempath)
    call mem_setptr(this%depth, 'DEPTH', this%input_mempath)
    !
    ! -- checkin EVT input context pointers
    call mem_checkin(this%surface, 'SURFACE', this%memoryPath, &
                     'SURFACE', this%input_mempath)
    call mem_checkin(this%rate, 'RATE', this%memoryPath, &
                     'RATE', this%input_mempath)
    call mem_checkin(this%depth, 'DEPTH', this%memoryPath, &
                     'DEPTH', this%input_mempath)
    !
    ! -- set list input segment descriptors
    if (.not. this%read_as_arrays) then
      if (this%nseg > 1) then
        !
        ! -- set pxdp and petm input context pointers
        call mem_setptr(this%pxdp, 'PXDP', this%input_mempath)
        call mem_setptr(this%petm, 'PETM', this%input_mempath)
        !
        ! -- checkin pxdp and petm input context pointers
        call mem_checkin(this%pxdp, 'PXDP', this%memoryPath, &
                         'PXDP', this%input_mempath)
        call mem_checkin(this%petm, 'PETM', this%memoryPath, &
                         'PETM', this%input_mempath)
      end if
      !
      if (this%surfratespecified) then
        !
        ! -- set petm0 input context pointer
        call mem_setptr(this%petm0, 'PETM0', this%input_mempath)
        !
        ! -- cehckin petm0 input context pointer
        call mem_checkin(this%petm0, 'PETM0', this%memoryPath, &
                         'PETM0', this%input_mempath)
      end if
    end if
  end subroutine evt_allocate_arrays

  !> @brief Source options specific to EvtType
  !<
  subroutine evt_source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(EvtType), intent(inout) :: this
    ! -- local
    logical(LGP) :: found_fixed_cell = .false.
    logical(LGP) :: found_readasarrays = .false.
    logical(LGP) :: found_surfratespec = .false.
    !
    ! -- source common bound options
    call this%BndExtType%source_options()
    !
    ! -- update defaults with idm sourced values
    call mem_set_value(this%fixed_cell, 'FIXED_CELL', &
                       this%input_mempath, found_fixed_cell)
    call mem_set_value(this%read_as_arrays, 'READASARRAYS', &
                       this%input_mempath, found_readasarrays)
    call mem_set_value(this%surfratespecified, 'SURFRATESPEC', &
                       this%input_mempath, found_surfratespec)
    !
    if (found_readasarrays) then
      if (this%dis%supports_layers()) then
        this%text = texta
      else
        errmsg = 'READASARRAYS option is not compatible with selected'// &
                 ' discretization type.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
    end if
    !
    if (found_readasarrays .and. found_surfratespec) then
      if (this%read_as_arrays) then
        errmsg = 'READASARRAYS option is not compatible with the'// &
                 ' SURF_RATE_SPECIFIED option.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
    end if
    !
    ! -- log evt specific options
    call this%evt_log_options(found_fixed_cell, found_readasarrays, &
                              found_surfratespec)
  end subroutine evt_source_options

  !> @brief Source options specific to EvtType
  !<
  subroutine evt_log_options(this, found_fixed_cell, found_readasarrays, &
                             found_surfratespec)
    ! -- modules
    use MemoryManagerModule, only: mem_reallocate, mem_setptr
    use MemoryManagerExtModule, only: mem_set_value
    use CharacterStringModule, only: CharacterStringType
    ! -- dummy
    class(EvtType), intent(inout) :: this
    logical(LGP), intent(in) :: found_fixed_cell
    logical(LGP), intent(in) :: found_readasarrays
    logical(LGP), intent(in) :: found_surfratespec
    ! -- formats
    character(len=*), parameter :: fmtihact = &
      &"(4x, 'EVAPOTRANSPIRATION WILL BE APPLIED TO HIGHEST ACTIVE CELL.')"
    character(len=*), parameter :: fmtfixedcell = &
      &"(4x, 'EVAPOTRANSPIRATION WILL BE APPLIED TO SPECIFIED CELL.')"
    character(len=*), parameter :: fmtreadasarrays = &
      &"(4x, 'EVAPOTRANSPIRATION INPUT WILL BE READ AS ARRAYS.')"
    character(len=*), parameter :: fmtsrz = &
      &"(4x, 'ET RATE AT SURFACE WILL BE ZERO.')"
    character(len=*), parameter :: fmtsrs = &
      &"(4x, 'ET RATE AT SURFACE WILL BE AS SPECIFIED BY PETM0.')"
    !
    ! -- log found options
    write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
      //' OPTIONS'
    !
    if (found_fixed_cell) then
      write (this%iout, fmtfixedcell)
    end if
    !
    if (found_readasarrays) then
      write (this%iout, fmtreadasarrays)
    end if
    !
    if (found_surfratespec) then
      write (this%iout, fmtsrs)
    end if
    !
    ! -- close logging block
    write (this%iout, '(1x,a)') &
      'END OF '//trim(adjustl(this%text))//' OPTIONS'
  end subroutine evt_log_options

  !> @brief Source the dimensions for this package
  !<
  subroutine evt_source_dimensions(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(EvtType), intent(inout) :: this
    ! -- local
    logical(LGP) :: found_nseg = .false.
    ! -- format
    character(len=*), parameter :: fmtnsegerr = &
      &"('Error: In EVT, NSEG must be > 0 but is specified as ',i0)"
    !
    ! Dimensions block is not required if:
    !   (1) discretization is DIS or DISV, and
    !   (2) READASARRAYS option has been specified.
    if (this%read_as_arrays) then
      this%maxbound = this%dis%get_ncpl()
      !
      ! -- verify dimensions were set
      if (this%maxbound <= 0) then
        write (errmsg, '(a)') &
          'MAXBOUND must be an integer greater than zero.'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
      !
    else
      !
      ! -- source maxbound
      call this%BndExtType%source_dimensions()
      !
      ! -- log found options
      write (this%iout, '(/1x,a)') 'PROCESSING '//trim(adjustl(this%text)) &
        //' DIMENSIONS'
      !
      ! -- update defaults with idm sourced values
      call mem_set_value(this%nseg, 'NSEG', this%input_mempath, found_nseg)
      !
      if (found_nseg) then
        !
        write (this%iout, '(4x,a,i0)') 'NSEG = ', this%nseg
        !
        if (this%nseg < 1) then
          write (errmsg, fmtnsegerr) this%nseg
          call store_error(errmsg)
          call store_error_filename(this%input_fname)
          !
        elseif (this%nseg > 1) then
          ! NSEG>1 is supported only if readasarrays is false
          if (this%read_as_arrays) then
            errmsg = 'In the EVT package, NSEG cannot be greater than 1'// &
                     ' when READASARRAYS is used.'
            call store_error(errmsg)
            call store_error_filename(this%input_fname)
          end if
          !
        end if
      end if
      !
      ! -- close logging block
      write (this%iout, '(1x,a)') &
        'END OF '//trim(adjustl(this%text))//' DIMENSIONS'
      !
    end if
    !
    ! -- Call define_listlabel to construct the list label that is written
    !    when PRINT_INPUT option is used.
    call this%define_listlabel()
  end subroutine evt_source_dimensions

  !> @brief Part of allocate and read
  !!
  !! If READASARRAYS has been specified, assign default IEVT = 1
  !<
  subroutine evt_read_initial_attr(this)
    ! -- dummy
    class(EvtType), intent(inout) :: this
    !
    if (this%read_as_arrays) then
      call this%default_nodelist()
    end if
  end subroutine evt_read_initial_attr

  !> @brief Read and Prepare
  !!
  !! Read itmp and new boundaries if itmp > 0
  !<
  subroutine evt_rp(this)
    use TdisModule, only: kper
    implicit none
    ! -- dummy
    class(EvtType), intent(inout) :: this
    !
    if (this%iper /= kper) return
    !
    if (this%read_as_arrays) then
      !
      ! -- update nodelist based on IEVT input
      call nodelist_update(this%nodelist, this%nbound, this%maxbound, &
                           this%dis, this%input_mempath)
      !
    else
      !
      ! -- process the input list arrays
      call this%BndExtType%bnd_rp()
      !
      ! -- ensure pxdp is monotonically increasing
      if (this%nseg > 1) then
        call this%check_pxdp()
      end if
      !
      ! -- Write the list to iout if requested
      if (this%iprpak /= 0) then
        call this%write_list()
      end if
      !
    end if
    !
    ! -- copy nodelist to nodesontop if not fixed cell
    if (.not. this%fixed_cell) call this%set_nodesontop()
  end subroutine evt_rp

  !> @brief Subroutine to check pxdp
  !!
  !! If the number of EVT segments (nseg) is greater than one, then
  !! pxdp must be monotically increasing from zero to one.  Check
  !! to make sure this is the case.
  !<
  subroutine check_pxdp(this)
    ! -- dummy
    class(EvtType), intent(inout) :: this !< EvtType
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: node
    integer(I4B) :: i
    integer(I4B) :: ierrmono
    real(DP) :: pxdp1
    real(DP) :: pxdp2
    character(len=15) :: nodestr
    ! -- formats
    character(len=*), parameter :: fmtpxdp0 = &
      &"('PXDP must be between 0 and 1.  Found ', G0, ' for cell ', A)"
    character(len=*), parameter :: fmtpxdp = &
      &"('PXDP is not monotonically increasing for cell ', A)"
    !
    ! -- check and make sure that pxdp is monotonically increasing and
    !    that pxdp values are between 0 and 1
    do n = 1, this%nbound
      node = this%nodelist(n)
      pxdp1 = DZERO
      ierrmono = 0
      segloop: do i = 1, this%nseg
        !
        ! -- set and check pxdp2
        if (i < this%nseg) then
          pxdp2 = this%pxdp(i, n)
          if (pxdp2 <= DZERO .or. pxdp2 >= DONE) then
            call this%dis%noder_to_string(node, nodestr)
            write (errmsg, fmtpxdp0) pxdp2, trim(nodestr)
            call store_error(errmsg)
          end if
        else
          pxdp2 = DONE
        end if
        !
        ! -- check for monotonically increasing condition
        if (pxdp2 - pxdp1 < DZERO) then
          if (ierrmono == 0) then
            ! -- only store mono error once for each node
            call this%dis%noder_to_string(node, nodestr)
            write (errmsg, fmtpxdp) trim(nodestr)
            call store_error(errmsg)
          end if
          ierrmono = 1
        end if
        pxdp1 = pxdp2
      end do segloop
    end do
    !
    ! -- terminate if errors encountered
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine check_pxdp

  !> @brief Store nodelist in nodesontop
  !<
  subroutine set_nodesontop(this)
    ! -- dummy
    class(EvtType), intent(inout) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- allocate if necessary
    if (.not. associated(this%nodesontop)) then
      allocate (this%nodesontop(this%maxbound))
    end if
    !
    ! -- copy nodelist into nodesontop
    do n = 1, this%nbound
      this%nodesontop(n) = this%nodelist(n)
    end do
  end subroutine set_nodesontop

  !> @brief Formulate the HCOF and RHS terms
  !<
  subroutine evt_cf(this)
    ! -- dummy
    class(EvtType) :: this
    ! -- local
    integer(I4B) :: i, iseg, node
    integer(I4B) :: idxdepth, idxrate
    real(DP) :: c, d, h, s, x
    real(DP) :: petm0
    real(DP) :: petm1, petm2, pxdp1, pxdp2, thcof, trhs
    !
    ! -- Return if no ET nodes
    if (this%nbound == 0) return
    !
    ! -- Calculate hcof and rhs for each ET node
    do i = 1, this%nbound
      !
      ! -- Find the node number
      if (this%fixed_cell) then
        node = this%nodelist(i)
      else
        node = this%nodesontop(i)
      end if
      !
      ! -- cycle if nonexistent bound
      if (node <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if
      !
      ! -- reset nodelist to highest active
      if (.not. this%fixed_cell) then
        if (this%ibound(node) == 0) &
          call this%dis%highest_active(node, this%ibound)
        this%nodelist(i) = node
      end if
      !
      ! -- set rhs and hcof to zero
      this%rhs(i) = DZERO
      this%hcof(i) = DZERO
      !
      ! -- if ibound is positive and not overlain by a lake, then add terms
      if (this%ibound(node) > 0 .and. this%ibound(node) /= IWETLAKE) then
        !
        if (this%iauxmultcol > 0) then
          c = this%rate(i) * this%dis%get_area(node) * &
              this%auxvar(this%iauxmultcol, i)
        else
          c = this%rate(i) * this%dis%get_area(node)
        end if
        s = this%surface(i)
        h = this%xnew(node)
        if (this%surfratespecified) then
          petm0 = this%petm0(i)
        end if
        !
        ! -- If head in cell is greater than or equal to SURFACE, ET is constant
        if (h >= s) then
          if (this%surfratespecified) then
            ! -- Subtract -PETM0 * max rate from RHS
            this%rhs(i) = this%rhs(i) + petm0 * c
          else
            ! -- Subtract -RATE from RHS
            this%rhs(i) = this%rhs(i) + c
          end if
        else
          ! -- If depth to water >= extinction depth, then ET is 0
          d = S - h
          x = this%depth(i)
          if (d < x) then
            ! -- Variable range. add ET terms to both RHS and HCOF.
            if (this%nseg > 1) then
              ! -- Determine which segment applies based on head, and
              !    calculate terms to add to RHS and HCOF
              !
              ! -- Set proportions corresponding to surface elevation
              !    and initial indices for depth and rate.
              ! -- Idxdepth will point to the elements of bound containing
              !    proportion of depth at the bottom of each segment.
              !    Idxrate will point to the elements of bound containing
              !    proportion of ET rate at the bottom of each segment.
              !    If surfratespecified is true, skip over the elements
              !    containing pxdp0 (=0.0) and petm0.
              pxdp1 = DZERO
              if (this%surfratespecified) then
                petm1 = petm0
              else
                petm1 = DONE
              end if
              ! -- Initialize indices to point to elements preceding
              !    pxdp1 and petm1 (values for lower end of segment 1).
              idxdepth = 0
              idxrate = 0
              ! -- Iterate through segments to find segment that contains
              !    current depth of head below ET surface.
              segloop: do iseg = 1, this%nseg
                ! -- Set proportions corresponding to lower end of
                !    segment
                if (iseg < this%nseg) then
                  ! -- Increment the indices for depth and rate
                  idxdepth = idxdepth + 1
                  idxrate = idxrate + 1
                  ! -- Get proportions for lower end of segment
                  pxdp2 = this%pxdp(idxdepth, i)
                  petm2 = this%petm(idxrate, i)
                else
                  pxdp2 = DONE
                  petm2 = DZERO
                end if
                if (d <= pxdp2 * x) then
                  ! -- head is in domain of this segment
                  exit segloop
                end if
                ! -- Proportions at lower end of segment will be for
                !    upper end of segment next time through loop
                pxdp1 = pxdp2
                petm1 = petm2
              end do segloop
              ! -- Calculate terms to add to RHS and HCOF based on
              !    segment that applies at head elevation
              thcof = -(petm1 - petm2) * c / ((pxdp2 - pxdp1) * x)
              trhs = thcof * (s - pxdp1 * x) + petm1 * c
            else
              ! -- Calculate terms to add to RHS and HCOF based on simple
              !    linear relation of ET vs. head for single segment
              trhs = c - c * s / x
              thcof = -c / x
            end if
            this%rhs(i) = this%rhs(i) + trhs
            this%hcof(i) = this%hcof(i) + thcof
          end if
        end if
      end if
      !
    end do
  end subroutine evt_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine evt_fc(this, rhs, ia, idxglo, matrix_sln)
    ! -- dummy
    class(EvtType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! -- local
    integer(I4B) :: i, n, ipos
    !
    ! -- Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      if (n <= 0) cycle
      ! -- reset hcof and rhs for excluded cells
      if (this%ibound(n) == IWETLAKE) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
    end do
  end subroutine evt_fc

  !> @brief Deallocate
  !<
  subroutine evt_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(EvtType) :: this
    !
    ! -- arrays
    if (associated(this%nodesontop)) deallocate (this%nodesontop)
    call mem_deallocate(this%surface, 'SURFACE', this%memoryPath)
    call mem_deallocate(this%rate, 'RATE', this%memoryPath)
    call mem_deallocate(this%depth, 'DEPTH', this%memoryPath)
    !
    if (.not. this%read_as_arrays) then
      if (this%nseg > 1) then
        call mem_deallocate(this%pxdp, 'PXDP', this%memoryPath)
        call mem_deallocate(this%petm, 'PETM', this%memoryPath)
      end if
      !
      if (this%surfratespecified) then
        call mem_deallocate(this%petm0, 'PETM0', this%memoryPath)
      end if
    end if
    !
    ! -- scalars
    call mem_deallocate(this%nseg)
    deallocate (this%segsdefined)
    deallocate (this%fixed_cell)
    deallocate (this%read_as_arrays)
    deallocate (this%surfratespecified)
    !
    ! -- Deallocate parent package
    call this%BndExtType%bnd_da()
  end subroutine evt_da

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !! option is used
  !<
  subroutine evt_define_listlabel(this)
    ! -- dummy
    class(EvtType), intent(inout) :: this
    ! -- local
    integer(I4B) :: nsegm1, i
    !
    ! -- create the header list label
    this%listlabel = trim(this%filtyp)//' NO.'
    if (this%dis%ndim == 3) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'ROW'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'COL'
    elseif (this%dis%ndim == 2) then
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'LAYER'
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'CELL2D'
    else
      write (this%listlabel, '(a, a7)') trim(this%listlabel), 'NODE'
    end if
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'SURFACE'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'MAX. RATE'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'EXT. DEPTH'
    !
    ! -- add headings for as many PXDP and PETM columns as needed
    nsegm1 = this%nseg - 1
    if (nsegm1 > 0) then
      do i = 1, nsegm1
        write (this%listlabel, '(a, a16)') trim(this%listlabel), 'PXDP'
      end do
      do i = 1, nsegm1
        write (this%listlabel, '(a, a16)') trim(this%listlabel), 'PETM'
      end do
    end if
    !
    ! -- PETM0, if SURF_RATE_SPECIFIED is used
    if (this%surfratespecified) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'PETM0'
    end if
    !
!    ! -- multiplier
!    if(this%multindex > 0) &
!      write(this%listlabel, '(a, a16)') trim(this%listlabel), 'MULTIPLIER'
    !
    ! -- boundary name
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
  end subroutine evt_define_listlabel

  !> @brief Assign default nodelist when READASARRAYS is specified.
  !!
  !! Equivalent to reading IEVT as CONSTANT 1
  !<
  subroutine default_nodelist(this)
    ! -- modules
    use SimModule, only: store_error
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(EvtType) :: this
    ! -- local
    integer(I4B) :: il, ir, ic, ncol, nrow, nlay, nodeu, noder, ipos
    !
    ! -- set variables
    if (this%dis%ndim == 3) then
      nlay = this%dis%mshape(1)
      nrow = this%dis%mshape(2)
      ncol = this%dis%mshape(3)
    elseif (this%dis%ndim == 2) then
      nlay = this%dis%mshape(1)
      nrow = 1
      ncol = this%dis%mshape(2)
    end if
    !
    ! -- Populate nodelist
    ipos = 1
    il = 1
    do ir = 1, nrow
      do ic = 1, ncol
        nodeu = get_node(il, ir, ic, nlay, nrow, ncol)
        noder = this%dis%get_nodenumber(nodeu, 0)
        this%nodelist(ipos) = noder
        ipos = ipos + 1
      end do
    end do
    !
    ! -- assign nbound.
    this%nbound = ipos - 1
    !
    ! -- if fixed_cell option not set, then need to store nodelist
    !    in the nodesontop array
    if (.not. this%fixed_cell) call this%set_nodesontop()
  end subroutine default_nodelist

  ! -- Procedures related to observations

  !> @brief Return true because EVT package supports observations
  !!
  !! Overrides BndType%bnd_obs_supported()
  !<
  logical function evt_obs_supported(this)
    ! -- dummy
    class(EvtType) :: this
    !
    evt_obs_supported = .true.
  end function evt_obs_supported

  !> @brief Store observation type supported by EVT package
  !!
  !! Overrides BndType%bnd_df_obs
  !<
  subroutine evt_df_obs(this)
    ! -- dummy
    class(EvtType) :: this
    ! -- local
    integer(I4B) :: indx
    !
    call this%obs%StoreObsType('evt', .true., indx)
    this%obs%obsData(indx)%ProcessIdPtr => DefaultObsIdProcessor
  end subroutine evt_df_obs

  !> @brief Return requested boundary value
  !<
  function evt_bound_value(this, col, row) result(bndval)
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy variables
    class(EvtType), intent(inout) :: this !< BndExtType object
    integer(I4B), intent(in) :: col
    integer(I4B), intent(in) :: row
    ! -- result
    real(DP) :: bndval
    ! -- local
    integer(I4B) :: idx
    !
    ! -- initialize
    idx = 0
    !
    select case (col)
    case (1)
      bndval = this%surface(row)
    case (2)
      if (this%iauxmultcol > 0) then
        bndval = this%rate(row) * this%auxvar(this%iauxmultcol, row)
      else
        bndval = this%rate(row)
      end if
    case (3)
      bndval = this%depth(row)
    case default
      if (col > 0) then
        if (this%nseg > 1) then
          if (col < (3 + this%nseg)) then
            idx = col - 3
            bndval = this%pxdp(idx, row)
          else if (col < (3 + (2 * this%nseg) - 1)) then
            idx = col - (3 + this%nseg - 1)
            bndval = this%petm(idx, row)
          else if (col == (3 + (2 * this%nseg) - 1)) then
            if (this%surfratespecified) then
              idx = 1
              bndval = this%petm0(row)
            end if
          end if
        else if (this%surfratespecified) then
          if (col == 4) then
            idx = 1
            bndval = this%petm0(row)
          end if
        end if
      end if
      !
      ! -- set error if idx not found
      if (idx == 0) then
        write (errmsg, '(a,i0,a)') &
          'Programming error. EVT bound value requested column '&
          &'outside range of ncolbnd (', this%ncolbnd, ').'
        call store_error(errmsg)
        call store_error_filename(this%input_fname)
      end if
      !
    end select
  end function evt_bound_value

  !> @brief Update the nodelist based on IEVT input
  !!
  !! This is a module scoped routine to check for IEVT input. If array input
  !! was provided, INIEVT and IEVT will be allocated in the input context.
  !! If the read state variable INIEVT is set to 1 during this period update,
  !! IEVT input was read and is used here to update the nodelist.
  !<
  subroutine nodelist_update(nodelist, nbound, maxbound, &
                             dis, input_mempath)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use BaseDisModule, only: DisBaseType
    ! -- dummy
    integer(I4B), dimension(:), contiguous, &
      pointer, intent(inout) :: nodelist
    class(DisBaseType), pointer, intent(in) :: dis
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(inout) :: nbound
    integer(I4B), intent(in) :: maxbound
    ! -- format
    character(len=24) :: aname = '     LAYER OR NODE INDEX'
    ! -- local
    integer(I4B), dimension(:), contiguous, pointer :: ievt => null()
    integer(I4B), pointer :: inievt => NULL()
    !
    ! -- set pointer to input context INIEVT
    call mem_setptr(inievt, 'INIEVT', input_mempath)
    !
    ! -- check INIEVT read state
    if (inievt == 1) then
      ! -- ievt was read this period
      !
      ! -- set pointer to input context IEVT
      call mem_setptr(ievt, 'IEVT', input_mempath)
      !
      ! -- update nodelist
      call dis%nlarray_to_nodelist(ievt, nodelist, &
                                   maxbound, nbound, aname)
    end if
  end subroutine nodelist_update

end module EvtModule

