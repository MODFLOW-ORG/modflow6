
module GwfHfbModule

  use KindModule, only: DP, I4B
  use SimVariablesModule, only: errmsg
  use SimModule, only: store_error, count_errors, store_error_filename
  use Xt3dModule, only: Xt3dType
  use GwfVscModule, only: GwfVscType
  use NumericalPackageModule, only: NumericalPackageType
  use BaseDisModule, only: DisBaseType
  use MatrixBaseModule

  implicit none

  private
  public :: GwfHfbType
  public :: hfb_cr

  type, extends(NumericalPackageType) :: GwfHfbType

    type(GwfVscType), pointer :: vsc => null() !< viscosity object
    integer(I4B), pointer :: maxhfb => null() !< max number of hfb's
    integer(I4B), pointer :: nhfb => null() !< number of hfb's
    integer(I4B), dimension(:), pointer, contiguous :: noden => null() !< first cell
    integer(I4B), dimension(:), pointer, contiguous :: nodem => null() !< second cell
    integer(I4B), dimension(:), pointer, contiguous :: idxloc => null() !< position in model ja
    real(DP), dimension(:), pointer, contiguous :: hydchr => null() !< hydraulic characteristic of the barrier
    real(DP), dimension(:), pointer, contiguous :: csatsav => null() !< value of condsat prior to hfb modification
    real(DP), dimension(:), pointer, contiguous :: condsav => null() !< saved conductance of combined npf and hfb
    type(Xt3dType), pointer :: xt3d => null() !< pointer to xt3d object
    !
    integer(I4B), dimension(:), pointer, contiguous :: ibound => null() !< pointer to model ibound
    integer(I4B), dimension(:), pointer, contiguous :: icelltype => null() !< pointer to model icelltype
    integer(I4B), dimension(:), pointer, contiguous :: ihc => null() !< pointer to model ihc
    integer(I4B), dimension(:), pointer, contiguous :: ia => null() !< pointer to model ia
    integer(I4B), dimension(:), pointer, contiguous :: ja => null() !< pointer to model ja
    integer(I4B), dimension(:), pointer, contiguous :: jas => null() !< pointer to model jas
    integer(I4B), dimension(:), pointer, contiguous :: isym => null() !< pointer to model isym
    real(DP), dimension(:), pointer, contiguous :: condsat => null() !< pointer to model condsat
    real(DP), dimension(:), pointer, contiguous :: top => null() !< pointer to model top
    real(DP), dimension(:), pointer, contiguous :: bot => null() !< pointer to model bot
    real(DP), dimension(:), pointer, contiguous :: hwva => null() !< pointer to model hwva
    real(DP), dimension(:), pointer, contiguous :: hnew => null() !< pointer to model xnew
    !
    ! -- viscosity flag
    integer(I4B), pointer :: ivsc => null() !< flag indicating if viscosity is active in the model

  contains

    procedure :: hfb_ar
    procedure :: hfb_rp
    procedure :: hfb_fc
    procedure :: hfb_cq
    procedure :: hfb_da
    procedure :: allocate_scalars
    procedure, private :: allocate_arrays
    procedure, private :: source_options
    procedure, private :: source_dimensions
    procedure, private :: source_data
    procedure, private :: check_data
    procedure, private :: condsat_reset
    procedure, private :: condsat_modify

  end type GwfHfbType

contains

  !> @brief Create a new hfb object
  !<
  subroutine hfb_cr(hfbobj, name_model, input_mempath, inunit, iout)
    ! -- dummy
    type(GwfHfbType), pointer :: hfbobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    !
    ! -- Create the object
    allocate (hfbobj)
    !
    ! -- create name and memory path
    call hfbobj%set_names(1, name_model, 'HFB', 'HFB', input_mempath)
    !
    ! -- Allocate scalars
    call hfbobj%allocate_scalars()
    !
    ! -- Save unit numbers
    hfbobj%inunit = inunit
    hfbobj%iout = iout
  end subroutine hfb_cr

  !> @brief Allocate and read
  !<
  subroutine hfb_ar(this, ibound, xt3d, dis, invsc, vsc)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use MemoryHelperModule, only: create_mem_path
    ! -- dummy
    class(GwfHfbType) :: this
    integer(I4B), dimension(:), pointer, contiguous :: ibound
    type(Xt3dType), pointer :: xt3d
    class(DisBaseType), pointer, intent(inout) :: dis !< discretization package
    integer(I4B), pointer :: invsc !< indicates if viscosity package is active
    type(GwfVscType), pointer, intent(in) :: vsc !< viscosity package
    ! -- formats
    character(len=*), parameter :: fmtheader = &
      "(1x, /1x, 'HFB -- HORIZONTAL FLOW BARRIER PACKAGE, VERSION 8, ', &
      &'4/24/2015 INPUT READ FROM MEMPATH: ', a, /)"
    !
    ! -- Print a message identifying the node property flow package.
    write (this%iout, fmtheader) this%input_mempath
    !
    ! -- Set pointers
    this%dis => dis
    this%ibound => ibound
    this%xt3d => xt3d
    !
    call mem_setptr(this%icelltype, 'ICELLTYPE', &
                    create_mem_path(this%name_model, 'NPF'))
    call mem_setptr(this%ihc, 'IHC', create_mem_path(this%name_model, 'CON'))
    call mem_setptr(this%ia, 'IA', create_mem_path(this%name_model, 'CON'))
    call mem_setptr(this%ja, 'JA', create_mem_path(this%name_model, 'CON'))
    call mem_setptr(this%jas, 'JAS', create_mem_path(this%name_model, 'CON'))
    call mem_setptr(this%isym, 'ISYM', create_mem_path(this%name_model, 'CON'))
    call mem_setptr(this%condsat, 'CONDSAT', create_mem_path(this%name_model, &
                                                             'NPF'))
    call mem_setptr(this%top, 'TOP', create_mem_path(this%name_model, 'DIS'))
    call mem_setptr(this%bot, 'BOT', create_mem_path(this%name_model, 'DIS'))
    call mem_setptr(this%hwva, 'HWVA', create_mem_path(this%name_model, 'CON'))
    !
    call this%source_options()
    call this%source_dimensions()
    call this%allocate_arrays()
    !
    ! --  If vsc package active, set ivsc
    if (invsc /= 0) then
      this%ivsc = 1
      this%vsc => vsc
      !
      ! -- Notify user via listing file viscosity accounted for by HFB package
      write (this%iout, '(/1x,a,a)') 'Viscosity active in ', &
        trim(this%filtyp)//' Package calculations: '//trim(adjustl(this%packName))
    end if
  end subroutine hfb_ar

  !> @brief Check for new HFB stress period data
  !<
  subroutine hfb_rp(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use TdisModule, only: kper
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    integer(I4B), pointer :: iper
    ! -- formats
    character(len=*), parameter :: fmtlsp = &
      &"(1X,/1X,'REUSING ',A,'S FROM LAST STRESS PERIOD')"

    call mem_setptr(iper, 'IPER', this%input_mempath)
    if (iper == kper) then
      call this%condsat_reset()
      call this%source_data()
      call this%condsat_modify()
    else
      write (this%iout, fmtlsp) 'HFB'
    end if
  end subroutine hfb_rp

  !> @brief Fill matrix terms
  !!
  !! Fill amatsln for the following conditions:
  !!   1. XT3D
  !!     OR
  !!   2. Not Newton, and
  !!   3. Cell type n is convertible or cell type m is convertible
  !<
  subroutine hfb_fc(this, kiter, matrix_sln, idxglo, rhs, hnew)
    ! -- modules
    use ConstantsModule, only: DHALF, DZERO, DONE
    ! -- dummy
    class(GwfHfbType) :: this
    integer(I4B) :: kiter
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(:) :: idxglo
    real(DP), intent(inout), dimension(:) :: rhs
    real(DP), intent(inout), dimension(:) :: hnew
    ! -- local
    integer(I4B) :: nodes, nja
    integer(I4B) :: ihfb, n, m
    integer(I4B) :: ipos
    integer(I4B) :: idiag, isymcon
    integer(I4B) :: ixt3d
    real(DP) :: cond, condhfb, aterm
    real(DP) :: fawidth, faheight
    real(DP) :: topn, topm, botn, botm
    real(DP) :: viscratio
    !
    ! -- initialize variables
    viscratio = DONE
    nodes = this%dis%nodes
    nja = this%dis%con%nja
    if (associated(this%xt3d%ixt3d)) then
      ixt3d = this%xt3d%ixt3d
    else
      ixt3d = 0
    end if
    !
    if (ixt3d > 0) then
      !
      do ihfb = 1, this%nhfb
        n = min(this%noden(ihfb), this%nodem(ihfb))
        m = max(this%noden(ihfb), this%nodem(ihfb))
        ! -- Skip if either cell is inactive.
        if (this%ibound(n) == 0 .or. this%ibound(m) == 0) cycle
        !!! if(this%icelltype(n) == 1 .or. this%icelltype(m) == 1) then
        if (this%ivsc /= 0) then
          call this%vsc%get_visc_ratio(n, m, hnew(n), hnew(m), viscratio)
        end if
        ! -- Compute scale factor for hfb correction
        if (this%hydchr(ihfb) > DZERO) then
          if (this%inewton == 0) then
            ipos = this%idxloc(ihfb)
            topn = this%top(n)
            topm = this%top(m)
            botn = this%bot(n)
            botm = this%bot(m)
            if (this%icelltype(n) == 1) then
              if (hnew(n) < topn) topn = hnew(n)
            end if
            if (this%icelltype(m) == 1) then
              if (hnew(m) < topm) topm = hnew(m)
            end if
            if (this%ihc(this%jas(ipos)) == 2) then
              faheight = min(topn, topm) - max(botn, botm)
            else
              faheight = DHALF * ((topn - botn) + (topm - botm))
            end if
            fawidth = this%hwva(this%jas(ipos))
            condhfb = this%hydchr(ihfb) * viscratio * &
                      fawidth * faheight
          else
            condhfb = this%hydchr(ihfb) * viscratio
          end if
        else
          condhfb = this%hydchr(ihfb)
        end if
        ! -- Make hfb corrections for xt3d
        call this%xt3d%xt3d_fhfb(kiter, nodes, nja, matrix_sln, idxglo, &
                                 rhs, hnew, n, m, condhfb)
      end do
      !
    else
      !
      ! -- For Newton, the effect of the barrier is included in condsat.
      if (this%inewton == 0) then
        do ihfb = 1, this%nhfb
          ipos = this%idxloc(ihfb)
          aterm = matrix_sln%get_value_pos(idxglo(ipos))
          n = this%noden(ihfb)
          m = this%nodem(ihfb)
          if (this%ibound(n) == 0 .or. this%ibound(m) == 0) cycle
          !
          if (this%ivsc /= 0) then
            call this%vsc%get_visc_ratio(n, m, hnew(n), hnew(m), viscratio)
          end if
          !
          if (this%icelltype(n) == 1 .or. this%icelltype(m) == 1 .or. &
              this%ivsc /= 0) then
            !
            ! -- Calculate hfb conductance
            topn = this%top(n)
            topm = this%top(m)
            botn = this%bot(n)
            botm = this%bot(m)
            if (this%icelltype(n) == 1) then
              if (hnew(n) < topn) topn = hnew(n)
            end if
            if (this%icelltype(m) == 1) then
              if (hnew(m) < topm) topm = hnew(m)
            end if
            if (this%ihc(this%jas(ipos)) == 2) then
              faheight = min(topn, topm) - max(botn, botm)
            else
              faheight = DHALF * ((topn - botn) + (topm - botm))
            end if
            if (this%hydchr(ihfb) > DZERO) then
              fawidth = this%hwva(this%jas(ipos))
              condhfb = this%hydchr(ihfb) * viscratio * &
                        fawidth * faheight
              cond = aterm * condhfb / (aterm + condhfb)
            else
              cond = -aterm * this%hydchr(ihfb)
            end if
            !
            ! -- Save cond for budget calculation
            this%condsav(ihfb) = cond
            !
            ! -- Fill row n diag and off diag
            idiag = this%ia(n)
            call matrix_sln%add_value_pos(idxglo(idiag), aterm - cond)
            call matrix_sln%set_value_pos(idxglo(ipos), cond)
            !
            ! -- Fill row m diag and off diag
            isymcon = this%isym(ipos)
            idiag = this%ia(m)
            call matrix_sln%add_value_pos(idxglo(idiag), aterm - cond)
            call matrix_sln%set_value_pos(idxglo(isymcon), cond)
            !
          end if
        end do
      end if
      !
    end if
  end subroutine hfb_fc

  !> @brief flowja will automatically include the effects of the hfb for
  !! confined and newton cases when xt3d is not used.
  !!
  !! This method recalculates flowja for the other cases.
  !<
  subroutine hfb_cq(this, hnew, flowja)
    ! -- modules
    use ConstantsModule, only: DHALF, DZERO, DONE
    ! -- dummy
    class(GwfHfbType) :: this
    real(DP), intent(inout), dimension(:) :: hnew
    real(DP), intent(inout), dimension(:) :: flowja
    ! -- local
    integer(I4B) :: ihfb, n, m
    integer(I4B) :: ipos
    real(DP) :: qnm
    real(DP) :: cond
    integer(I4B) :: ixt3d
    real(DP) :: condhfb
    real(DP) :: fawidth, faheight
    real(DP) :: topn, topm, botn, botm
    real(DP) :: viscratio
    !
    ! -- initialize viscratio
    viscratio = DONE
    !
    if (associated(this%xt3d%ixt3d)) then
      ixt3d = this%xt3d%ixt3d
    else
      ixt3d = 0
    end if
    !
    if (ixt3d > 0) then
      !
      do ihfb = 1, this%nhfb
        n = min(this%noden(ihfb), this%nodem(ihfb))
        m = max(this%noden(ihfb), this%nodem(ihfb))
        ! -- Skip if either cell is inactive.
        if (this%ibound(n) == 0 .or. this%ibound(m) == 0) cycle
        !!! if(this%icelltype(n) == 1 .or. this%icelltype(m) == 1) then
        if (this%ivsc /= 0) then
          call this%vsc%get_visc_ratio(n, m, hnew(n), hnew(m), viscratio)
        end if
        !
        ! -- Compute scale factor for hfb correction
        if (this%hydchr(ihfb) > DZERO) then
          if (this%inewton == 0) then
            ipos = this%idxloc(ihfb)
            topn = this%top(n)
            topm = this%top(m)
            botn = this%bot(n)
            botm = this%bot(m)
            if (this%icelltype(n) == 1) then
              if (hnew(n) < topn) topn = hnew(n)
            end if
            if (this%icelltype(m) == 1) then
              if (hnew(m) < topm) topm = hnew(m)
            end if
            if (this%ihc(this%jas(ipos)) == 2) then
              faheight = min(topn, topm) - max(botn, botm)
            else
              faheight = DHALF * ((topn - botn) + (topm - botm))
            end if
            fawidth = this%hwva(this%jas(ipos))
            condhfb = this%hydchr(ihfb) * viscratio * &
                      fawidth * faheight
          else
            condhfb = this%hydchr(ihfb)
          end if
        else
          condhfb = this%hydchr(ihfb)
        end if
        ! -- Make hfb corrections for xt3d
        call this%xt3d%xt3d_flowjahfb(n, m, hnew, flowja, condhfb)
      end do
      !
    else
      !
      ! -- Recalculate flowja for non-newton unconfined.
      if (this%inewton == 0) then
        do ihfb = 1, this%nhfb
          n = this%noden(ihfb)
          m = this%nodem(ihfb)
          if (this%ibound(n) == 0 .or. this%ibound(m) == 0) cycle
          if (this%icelltype(n) == 1 .or. this%icelltype(m) == 1 .or. &
              this%ivsc /= 0) then
            ipos = this%dis%con%getjaindex(n, m)
            !
            ! -- condsav already accnts for visc adjustment
            cond = this%condsav(ihfb)
            qnm = cond * (hnew(m) - hnew(n))
            flowja(ipos) = qnm
            ipos = this%dis%con%getjaindex(m, n)
            flowja(ipos) = -qnm
            !
          end if
        end do
      end if
      !
    end if
  end subroutine hfb_cq

  !> @brief Deallocate memory
  !<
  subroutine hfb_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(GwfHfbType) :: this
    !
    ! -- Scalars
    call mem_deallocate(this%maxhfb)
    call mem_deallocate(this%nhfb)
    call mem_deallocate(this%ivsc)
    !
    ! -- Arrays
    if (this%inunit > 0) then
      call mem_deallocate(this%noden)
      call mem_deallocate(this%nodem)
      call mem_deallocate(this%hydchr)
      call mem_deallocate(this%idxloc)
      call mem_deallocate(this%csatsav)
      call mem_deallocate(this%condsav)
    end if
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
    !
    ! -- nullify pointers
    this%xt3d => null()
    this%inewton => null()
    this%ibound => null()
    this%icelltype => null()
    this%ihc => null()
    this%ia => null()
    this%ja => null()
    this%jas => null()
    this%isym => null()
    this%condsat => null()
    this%top => null()
    this%bot => null()
    this%hwva => null()
    this%vsc => null()
  end subroutine hfb_da

  !> @brief Allocate package scalars
  !<
  subroutine allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfHfbType) :: this
    !
    ! -- allocate scalars in NumericalPackageType
    call this%NumericalPackageType%allocate_scalars()
    !
    ! -- allocate scalars
    call mem_allocate(this%maxhfb, 'MAXHFB', this%memoryPath)
    call mem_allocate(this%nhfb, 'NHFB', this%memoryPath)
    !
    ! -- allocate flag for determining if vsc active
    call mem_allocate(this%ivsc, 'IVSC', this%memoryPath)
    !
    ! -- initialize
    this%maxhfb = 0
    this%nhfb = 0
    this%ivsc = 0
  end subroutine allocate_scalars

  !> @brief Allocate package arrays
  !<
  subroutine allocate_arrays(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    integer(I4B) :: ihfb
    !
    call mem_allocate(this%noden, this%maxhfb, 'NODEN', this%memoryPath)
    call mem_allocate(this%nodem, this%maxhfb, 'NODEM', this%memoryPath)
    call mem_allocate(this%hydchr, this%maxhfb, 'HYDCHR', this%memoryPath)
    call mem_allocate(this%idxloc, this%maxhfb, 'IDXLOC', this%memoryPath)
    call mem_allocate(this%csatsav, this%maxhfb, 'CSATSAV', this%memoryPath)
    call mem_allocate(this%condsav, this%maxhfb, 'CONDSAV', this%memoryPath)
    !
    ! -- initialize idxloc to 0
    do ihfb = 1, this%maxhfb
      this%idxloc(ihfb) = 0
    end do
  end subroutine allocate_arrays

  !> @ brief Source hfb input options
  !<
  subroutine source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfHfbInputModule, only: GwfHfbParamFoundType
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    type(GwfHfbParamFoundType) :: found

    ! update options from input context
    call mem_set_value(this%iprpak, 'PRINT_INPUT', this%input_mempath, &
                       found%print_input)

    ! log options
    write (this%iout, '(1x,a)') 'PROCESSING HFB OPTIONS'
    if (found%print_input) then
      write (this%iout, '(4x,a)') &
        'THE LIST OF HFBS WILL BE PRINTED.'
    end if
    write (this%iout, '(1x,a)') 'END OF HFB OPTIONS'
  end subroutine source_options

  !> @ brief Source hfb input options
  !<
  subroutine source_dimensions(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    use GwfHfbInputModule, only: GwfHfbParamFoundType
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    type(GwfHfbParamFoundType) :: found

    ! update dimensions from input context
    call mem_set_value(this%maxhfb, 'MAXBOUND', this%input_mempath, &
                       found%maxbound)

    ! log dimensions
    write (this%iout, '(/1x,a)') 'PROCESSING HFB DIMENSIONS'
    write (this%iout, '(4x,a,i7)') 'MAXHFB = ', this%maxhfb
    write (this%iout, '(1x,a)') 'END OF HFB DIMENSIONS'

    ! check dimensions
    if (this%maxhfb <= 0) then
      write (errmsg, '(a)') &
        'MAXHFB must be specified with value greater than zero.'
      call store_error(errmsg)
      call store_error_filename(this%input_mempath)
    end if
  end subroutine source_dimensions

  !> @ brief source hfb period data
  !<
  subroutine source_data(this)
    ! -- modules
    use TdisModule, only: kper
    use ConstantsModule, only: LINELENGTH
    use MemoryManagerModule, only: mem_setptr
    use GeomUtilModule, only: get_node
    ! -- dummy
    class(GwfHfbType), intent(inout) :: this
    ! -- local
    integer(I4B), dimension(:, :), pointer, contiguous :: cellids1, cellids2
    integer(I4B), dimension(:), pointer :: cellid1, cellid2
    real(DP), dimension(:), pointer, contiguous :: hydchr
    character(len=LINELENGTH) :: nodenstr, nodemstr
    integer(I4B), pointer :: nbound
    integer(I4B) :: n, nodeu1, nodeu2, noder1, noder2
    ! -- formats
    character(len=*), parameter :: fmthfb = "(i10, 2a10, 1(1pg15.6))"

    ! set input context pointers
    call mem_setptr(nbound, 'NBOUND', this%input_mempath)
    call mem_setptr(cellids1, 'CELLID1', this%input_mempath)
    call mem_setptr(cellids2, 'CELLID2', this%input_mempath)
    call mem_setptr(hydchr, 'HYDCHR', this%input_mempath)

    ! set nhfb
    this%nhfb = nbound

    ! log data
    write (this%iout, '(//,1x,a)') 'READING HFB DATA'
    if (this%iprpak > 0) then
      write (this%iout, '(3a10, 1a15)') 'HFB NUM', 'CELL1', 'CELL2', &
        'HYDCHR'
    end if

    ! update state
    do n = 1, this%nhfb

      ! set cellid
      cellid1 => cellids1(:, n)
      cellid2 => cellids2(:, n)

      ! set node user
      if (this%dis%ndim == 1) then
        nodeu1 = cellid1(1)
        nodeu2 = cellid2(1)
      elseif (this%dis%ndim == 2) then
        nodeu1 = get_node(cellid1(1), 1, cellid1(2), &
                          this%dis%mshape(1), 1, &
                          this%dis%mshape(2))
        nodeu2 = get_node(cellid2(1), 1, cellid2(2), &
                          this%dis%mshape(1), 1, &
                          this%dis%mshape(2))
      else
        nodeu1 = get_node(cellid1(1), cellid1(2), cellid1(3), &
                          this%dis%mshape(1), &
                          this%dis%mshape(2), &
                          this%dis%mshape(3))
        nodeu2 = get_node(cellid2(1), cellid2(2), cellid2(3), &
                          this%dis%mshape(1), &
                          this%dis%mshape(2), &
                          this%dis%mshape(3))
      end if

      ! set nodes
      noder1 = this%dis%get_nodenumber(nodeu1, 1)
      noder2 = this%dis%get_nodenumber(nodeu2, 1)
      if (noder1 <= 0 .or. &
          noder2 <= 0) then
        cycle
      else
        this%noden(n) = noder1
        this%nodem(n) = noder2
      end if

      this%hydchr(n) = hydchr(n)

      ! print input if requested
      if (this%iprpak /= 0) then
        call this%dis%noder_to_string(this%noden(n), nodenstr)
        call this%dis%noder_to_string(this%nodem(n), nodemstr)
        write (this%iout, fmthfb) n, trim(adjustl(nodenstr)), &
          trim(adjustl(nodemstr)), this%hydchr(n)
      end if
    end do

    ! check errors
    if (count_errors() > 0) then
      call store_error('Errors encountered in HFB input file.')
      call store_error_filename(this%input_fname)
    end if

    ! finalize logging
    write (this%iout, '(3x,i0,a,i0)') this%nhfb, &
      ' HFBs READ FOR STRESS PERIOD ', kper
    write (this%iout, '(1x,a)') 'END READING HFB DATA'

    ! input data check
    call this%check_data()
  end subroutine source_data

  !> @brief Check for hfb's between two unconnected cells and write a warning
  !!
  !! Store ipos in idxloc
  !<
  subroutine check_data(this)
    ! -- modules
    use ConstantsModule, only: LINELENGTH
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    integer(I4B) :: ihfb, n, m
    integer(I4B) :: ipos
    character(len=LINELENGTH) :: nodenstr, nodemstr
    logical :: found
    ! -- formats
    character(len=*), parameter :: fmterr = "(1x, 'HFB no. ',i0, &
      &' is between two unconnected cells: ', a, ' and ', a)"
    character(len=*), parameter :: fmtverr = "(1x, 'HFB no. ',i0, &
      &' is between two cells not horizontally connected: ', a, ' and ', a)"
    !
    do ihfb = 1, this%nhfb
      n = this%noden(ihfb)
      m = this%nodem(ihfb)
      found = .false.
      do ipos = this%ia(n) + 1, this%ia(n + 1) - 1
        if (m == this%ja(ipos)) then
          found = .true.
          this%idxloc(ihfb) = ipos
          exit
        end if
      end do
      !
      ! -- check to make sure cells are connected
      if (.not. found) then
        call this%dis%noder_to_string(n, nodenstr)
        call this%dis%noder_to_string(m, nodemstr)
        write (errmsg, fmterr) ihfb, trim(adjustl(nodenstr)), &
          trim(adjustl(nodemstr))
        call store_error(errmsg)
      else
        !
        ! -- check to make sure cells are not vertically connected
        ipos = this%idxloc(ihfb)
        if (this%ihc(this%jas(ipos)) == 0) then
          call this%dis%noder_to_string(n, nodenstr)
          call this%dis%noder_to_string(m, nodemstr)
          write (errmsg, fmtverr) ihfb, trim(adjustl(nodenstr)), &
            trim(adjustl(nodemstr))
          call store_error(errmsg)
        end if
      end if
    end do
    !
    ! -- Stop if errors detected
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine check_data

  !> @brief Reset condsat to its value prior to being modified by hfb's
  !<
  subroutine condsat_reset(this)
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    integer(I4B) :: ihfb
    integer(I4B) :: ipos
    !
    do ihfb = 1, this%nhfb
      ipos = this%idxloc(ihfb)
      this%condsat(this%jas(ipos)) = this%csatsav(ihfb)
    end do
  end subroutine condsat_reset

  !> @brief Modify condsat
  !!
  !! Modify condsat for the following conditions:
  !!   1.  If Newton is active
  !!   2.  If icelltype for n and icelltype for m is 0
  !<
  subroutine condsat_modify(this)
    ! -- modules
    use ConstantsModule, only: DHALF, DZERO
    ! -- dummy
    class(GwfHfbType) :: this
    ! -- local
    integer(I4B) :: ihfb, n, m
    integer(I4B) :: ipos
    real(DP) :: cond, condhfb
    real(DP) :: fawidth, faheight
    real(DP) :: topn, topm, botn, botm
    !
    do ihfb = 1, this%nhfb
      ipos = this%idxloc(ihfb)
      cond = this%condsat(this%jas(ipos))
      this%csatsav(ihfb) = cond
      n = this%noden(ihfb)
      m = this%nodem(ihfb)
      !
      if (this%inewton == 1 .or. &
          (this%icelltype(n) == 0 .and. this%icelltype(m) == 0)) then
        !
        ! -- Calculate hfb conductance
        topn = this%top(n)
        topm = this%top(m)
        botn = this%bot(n)
        botm = this%bot(m)
        if (this%ihc(this%jas(ipos)) == 2) then
          faheight = min(topn, topm) - max(botn, botm)
        else
          faheight = DHALF * ((topn - botn) + (topm - botm))
        end if
        if (this%hydchr(ihfb) > DZERO) then
          fawidth = this%hwva(this%jas(ipos))
          condhfb = this%hydchr(ihfb) * &
                    fawidth * faheight
          cond = cond * condhfb / (cond + condhfb)
        else
          cond = -cond * this%hydchr(ihfb)
        end if
        this%condsat(this%jas(ipos)) = cond
      end if
    end do
  end subroutine condsat_modify

end module GwfHfbModule
