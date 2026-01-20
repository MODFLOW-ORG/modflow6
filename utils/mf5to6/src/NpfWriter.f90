module NpfWriterModule

  use ConstantsModule, only: DONE, DZERO, DEM7, DHALF, DTWO
  use ConstantsPHMFModule, only: HDRYDEFAULT, HNOFLODEFAULT
  use FileListModule, only: FileListType
  use FileTypeModule, only: FileType
  use FileWriterModule, only: FileWriterType
  use GLOBAL, only: ncol, nrow, nlay, ITMUNI, LENUNI, LAYCBD, NCNFBD, BOTM, &
                    LBOTM, cbcfilename
  use GlobalVariablesModule, only: echo
  use GWFBCFMODULE, only: LAYCON, WetdryBcf => WETDRY, vcont, TRPY
  use GWFLPFMODULE, only: LAYTYP, LAYWET, VKCB, WetdryLpf => WETDRY, HANI
  use GWFUPWMODULE, only: LaytypUpw, LaywetUpw => LAYWET, VkcbUpw => VKCB, &
                          WetdryUpw => WETDRY, HaniUpw => HANI
  use InputOutputModule, only: GetUnit, openfile
  use SimPHMFModule, only: count_errors, store_error, store_note, &
                           store_warning, ustop
  use UtilitiesModule, only: Write2dValues

  type, extends(FileWriterType) :: NpfWriterType
    double precision            :: Hnoflo = hnoflodefault
    double precision            :: Hdry   = hdrydefault
    double precision            :: Wetfct = DZERO
    double precision            :: MinSatThickness = DEM7
    double precision, allocatable, dimension(:,:,:) :: hk
    double precision, allocatable, dimension(:,:,:) :: vk
    double precision, allocatable, dimension(:,:,:) :: hani
    double precision, pointer,     dimension(:,:,:) :: WetDry => null()
    integer                     :: Inpfcb = 0
    integer                     :: Iwetit = 1
    integer                     :: Ihdwet = 0
    integer                     :: NumConvertible = 0
    integer, pointer            :: Nlaynew => null()
    integer, allocatable, dimension(:) :: Icelltype ! dim = nlay
    integer, pointer,     dimension(:) :: Layptr => null()
    logical                     :: VertFlowCorr = .false.
    logical                     :: ConstantCv = .false.
    logical                     :: CvCorrection = .false.
    logical                     :: Rewet = .false.
    logical                     :: Newton = .false.
    logical                     :: ThickStrt = .false.
    logical                     :: BottomHeadDamping = .false.
    logical                     :: Perched = .false.
    logical                     :: VariableCV = .false.
    logical                     :: Dewatered = .false.
    logical                     :: UseHani = .false.
    character(len=3)            :: FlowPackage = ''
    character(len=11)           :: CellAveraging = 'HARMONIC'
    type(FileListType), pointer :: Mf2005Files => null()
    type(FileListType), pointer :: Mf6Files => null()
  contains
    procedure :: InitializeFile
    procedure :: AllocateArrays
    procedure :: WriteFile
    procedure, private :: WriteOptions
    procedure, private :: assign_vk
  end type NpfWriterType

contains

  subroutine AllocateArrays(this)
    implicit none
    ! dummy
    class(NpfWriterType), intent(inout) :: this
    !
    if (.not. allocated(this%Icelltype)) then
      allocate(this%Icelltype(nlay))
      allocate(this%hk(ncol,nrow,nlay))
      allocate(this%vk(ncol,nrow,this%Nlaynew))
    endif
    this%hk = DZERO
    this%vk = DONE
    !
    return
  end subroutine AllocateArrays

  subroutine InitializeFile(this, fname, ftype, pkgname)
    implicit none
    ! dummy
    class(NpfWriterType), intent(inout) :: this
    character(len=*), intent(in)    :: fname
    character(len=*), intent(in)    :: ftype
    character(len=*), intent(in), optional :: pkgname
    ! local
    !
    if (present(pkgname)) then
      call this%FileWriterType%InitializeFile(fname, ftype, pkgname)
    else
      call this%FileWriterType%InitializeFile(fname, ftype)
    endif
    !
    return
  end subroutine InitializeFile

  subroutine WriteFile(this)
    implicit none
    ! dummy
    class(NpfWriterType), intent(inout) :: this
    ! local
    integer :: i, iprni, iprnr, iu, iwetdry, j, k, knew, &
               idim1, idim2, idim3
    character(len=500) :: msg
    logical :: constant, vkerr, wettable, writek33
    double precision :: hk, val0, valt
    ! formats
    5  format()
    30 format(4x,a,2x,i0)
    40 format(4x,a,2x,g16.9)
    50 format(2x,a)
    60 format(a)
    100 format(a,1x,i0,1x,a,1x,i0,1x,a,1x,g10.3)
    !
!    ! Rewetting settings
!    if (associated(WetdryBcf)) then
!      this%WetDry => WetdryBcf
!    elseif (associated(WetdryLpf)) then
!      this%WetDry => WetdryLpf
!    elseif (associated(WetdryUpw)) then
!      this%WetDry => WetdryUpw
!    endif
    if (this%FlowPackage == 'UPW') then
      this%Rewet = .false.
      this%WetDry = 0.01d0
    endif
    !
    if (echo) then
      iprni = 3
      iprnr = 12
    else
      iprni = -3
      iprnr = -12
    endif
    call this%WriteOptions()
    !
    iu = this%fileobj%IUnit
    write(iu,5)
    write(iu,60)'BEGIN GRIDDATA'
    !
    ! Icelltype (0: confined; >0: convertible; <0: see input instruct.)
    !           (will be constant for each layer)
    write(iu,50)'Icelltype LAYERED'
    do k=1,nlay
      write(iu,30)'CONSTANT', this%Icelltype(k)
      if (LAYCBD(k) /= 0) then
        ! Assume quasi-3d unit is not convertible
        write(iu,30)'CONSTANT', 0
      endif
    enddo
    !
    ! HK
    write(iu,50)'K  LAYERED'
    knew = 0
    do k=1,nlay
      knew = knew + 1
      constant = .true.
      val0 = this%hk(1,1,k)
      do i=1,nrow
        do j=1,ncol
          if (this%hk(j,i,k) /= val0) then
            constant = .false.
            exit
          endif
        enddo
      enddo
      if (constant) then
        write(iu,40)'CONSTANT', val0
      else
        write(iu,30)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
        call Write2dValues(iu,NROW,NCOL,this%hk(:,:,k))
      endif
      if (LAYCBD(k) /= 0) then
        knew = knew + 1
        ! Get a default value for HK
        hk = get_q3d_hk()
        write(iu,40)'CONSTANT', hk
        write(msg,100)'Original quasi-3D unit below layer',k,'is now layer', &
                      knew,'and has been assigned HK value:',hk
        call store_note(msg)
      endif
    enddo
    !
    ! VK (may be problematic when converting from BCF? May need to solve
    !     a system of NLAY-1 equations? What about Quasi-3D layers?)
    ! -- Don't support BCF when NLAY > 1.  For LPF, convert Q-3D layers to
    !    active cells.
    ! -- Better: For BCF, assign VK = HK for original BCF layers, and assign
    !    VK = (default HK)/100 for new layer replacing a quasi-3d unit
    if (this%FlowPackage == 'BCF') then
      ! Assign all VK values
      vkerr = .false.
      call this%assign_vk(vkerr)
      if (vkerr) then
        do k=1,nlay
          do i=1,nrow
            do j=1,ncol
              this%vk(j,i,k) = this%hk(j,i,k)
            enddo
          enddo
        enddo
        if (nlay>1) then
          msg = 'For conversion of BCF layers to NPF package, definition of' // &
                ' vertical hydraulic conductivity (VK) of each layer does' // &
                ' not have a unique solution. For each layer, VK has been' // &
                ' set to horizontal hydraulic conductivity (HK) of the' // &
                ' original layer. These values likely will need to be' // &
                ' adjusted by the modeler to approximate the original model.'
          call store_note(msg)
          if (NCNFBD > 0) then
            msg = 'For conversion of BCF quasi-3D unit(s) to NPF package,' // &
                  ' default values of HK and VK have been assigned. These' // &
                  ' values likely will need to be adjusted by the modeler' // &
                  ' to approximate the original model.'
            call store_note(msg)
          endif
        endif
      endif
    elseif (this%FlowPackage == 'LPF') then
      ! Ned todo: is anything needed here?
      ! Probably need to populate vk?
    elseif (this%FlowPackage == 'UPW') then
      ! Anything here?
    endif
    !
    ! -- Determine if necessary to write K33
    writek33 = .false.
    if (nlay > 1) then
      layerloop: do k=1,this%Nlaynew
        do i=1,nrow
          do j=1,ncol
            if(this%vk(j,i,k) /= this%hk(j,i,k)) then
              writek33 = .true.
              exit layerloop
            endif
          enddo
        enddo
      enddo layerloop
    endif
    !
    if (writek33) then    
      write(iu,50)'K33  LAYERED'
      do k=1,this%Nlaynew
        constant = .true.
        val0 = this%vk(1,1,k)
        do i=1,nrow
          do j=1,ncol
            if (this%vk(j,i,k) /= val0) then
              constant = .false.
              exit
            endif
          enddo
        enddo
        if (constant) then
          write(iu,40)'CONSTANT', val0
        else
          write(iu,30)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
          call Write2dValues(iu,NROW,NCOL,this%vk(:,:,k))
        endif
      enddo
    endif
    !
    ! WETDRY
    if (this%Rewet) then
      write(iu,50)'WETDRY LAYERED'
      iwetdry = 0 ! index to WETDRY array
      do k=1,nlay
        ! Determine if this layer is wettable
        wettable = .false.
        select case (this%FlowPackage)
        case ('LPF')
          if (LAYTYP(k) /= 0 .and. LAYWET(k) /= 0) then
            wettable = .true.
          endif
        case ('UPW')
          if (LAYTYPUPW(k) /= 0) then ! .and. LaywetUpw(k) /= 0) then
            wettable = .true.
          endif
        case ('BCF')
          if (LAYCON(k) == 1 .or. LAYCON(k) == 3) then
            wettable = .true.
          endif
        end select
        !
        ! If this layer is wettable, need nonzero Wetdry values
        if (wettable) then
          ! determine if laywet for this layer is constant
          constant = .true.
          val0 = this%WetDry(1,1,1)
          idim1 = size(this%WetDry,1)
          idim2 = size(this%WetDry,2)
          idim3 = size(this%WetDry,3)
          if (idim1 == ncol .and. idim2 == nrow .and. idim3 == nlay) then
            val0 = this%WetDry(1,1,k)
            do i=1,nrow
              do j=1,ncol
                if (this%WetDry(j,i,k) /= val0) then
                  constant = .false.
                  exit
                endif
              enddo
            enddo
          endif
          if (constant) then
            write(iu,40)'CONSTANT', val0
          else
            write(iu,30)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
            call Write2dValues(iu,NROW,NCOL,this%WetDry(:,:,k))
          endif
        else
          ! Layer is not wettable in LPF, but needs to be assigned Wetdry = 0 in NPF
          write(iu,40)'CONSTANT', DZERO
        endif
        if (LAYCBD(k) /= 0) then
          ! Need to write WETDRY array for quasi-3d unit.
          ! Quasi-3d confining bed would not be rewettable, so use Wetdry = 0
          write(iu,40)'CONSTANT', DZERO
        endif
      enddo
    endif
    !
    ! HANI -- convert to K22
    if (this%UseHani) then
      write(iu,50)'K22  LAYERED'
      knew = 0
      do k=1,nlay
        knew = knew + 1
        constant = .true.
        val0 = this%hani(1,1,k) * this%hk(1,1,k)
        do i=1,nrow
          do j=1,ncol
            valt = this%hani(j,i,k) * this%hk(j,i,k)
            if (valt /= val0) then
              constant = .false.
              exit
            endif
          enddo
        enddo
        if (constant) then
          write(iu,40)'CONSTANT', val0
        else
          ! calculate k22
          do i=1,nrow
            do j=1,ncol
              this%hani(j,i,k) = this%hani(j,i,k) * this%hk(j,i,k)
            end do
          end do
          write(iu,30)'INTERNAL  FACTOR  1.0  IPRN ',iprnr
          call Write2dValues(iu,NROW,NCOL,this%hani(:,:,k))
        endif
        if (LAYCBD(k) /= 0) then
          knew = knew + 1
          ! Assign a default value for HANI
          write(iu,40)'CONSTANT', 1.0d0
          write(msg,100)'Original quasi-3D unit below layer',k,'is now layer', &
                        knew,'and has been assigned HANI value:',1.0d0
          call store_note(msg)
        endif
      enddo
    endif
    !
    ! ANGLEX (presumably not needed)
    !
    write(iu,60)'END GRIDDATA'
    !
    return
  end subroutine WriteFile

  subroutine WriteOptions(this)
    implicit none
    ! dummy
    class(NpfWriterType), intent(inout) :: this
    ! local
    integer :: icbc, iu
    type(FileType), pointer :: cbcfil
    ! formats
    5  format()
    10 format('BEGIN Options')
    20 format(2x,a,2x,a)
    30 format(2x,a,2x,i0)
    40 format(2x,a,2x,g16.9)
    50 format(2x,a)
    60 format(a)
    70 format(2x,a,2x,a,2x,g16.9,2x,a,2x,i0,2x,a,2x,i0)   
    100 format('END Options')
    !
    if (this%Newton) then
      this%Perched = .false.
      this%VariableCV = .false.
      this%Dewatered = .false.
      this%BottomHeadDamping = .true.
    endif
    !
    iu = this%fileobj%IUnit
    ! Write BEGIN Options
    write(iu,5)
    write(iu,10)
    !
    ! Options related to vertical flow correction
    if (this%Perched) then
      write(iu,50)'PERCHED'
    endif
    if (this%VariableCV) then
      if (this%Dewatered) then
        write(iu,20)'VARIABLECV', 'DEWATERED'
      else
        write(iu,50)'VARIABLECV'
      endif
    endif
    !
    ! SAVE_FLOWS
    cbcfil => null()
    if (this%Inpfcb/=0) then
      icbc = this%Inpfcb
      if (icbc > 0) then
        cbcfil => this%Mf6Files%GetFileByUnit(icbc)
        if (associated(cbcfil)) then
          if (cbcfilename == '') then
            cbcfilename = cbcfil%FName
          endif
          if (cbcfilename .ne. '') then
            write(iu,20)'SAVE_FLOWS', trim(cbcfilename)
          endif
!          write(iu,50)'SAVE_FLOWS'
        endif
      endif
    endif
    !!
    !! HNOFLO
    !if (this%Hnoflo/=hnoflodefault) then
    !  write(line,40)'HNOFLO', this%Hnoflo
    !  write(iu,60) trim(line)
    !endif
    !!
    !! HDRY
    !if (this%Hdry/=hdrydefault) then
    !  write(line,40)'HDRY', this%Hdry
    !  write(iu,60) trim(line)
    !endif
    !
    ! CELL_AVERAGING
    if (this%CellAveraging /= 'HARMONIC') then
      write(iu,20)'ALTERNATIVE_CELL_AVERAGING', trim(this%CellAveraging)
    endif
    !
    ! CONSTANTCV
    if (this%ConstantCv) write(iu,50)'CONSTANTCV'
    !
    ! CVCORRECTION
    if (this%CvCorrection) write(iu,50)'CVCORRECTION'
    !
    ! REWET options
    if (this%Rewet) then
      write(iu,70) 'REWET', 'WETFCT', this%Wetfct, 'IWETIT', this%Iwetit,       &
                   'IHDWET', this%Ihdwet
      !write(iu,50)'REWET'
      !if (this%Wetfct/=DZERO) then
      !  write(line,40)'WETFCT',this%Wetfct
      !  write(iu,60) trim(line)
      !endif
      !if (this%Iwetit>1) then
      !  write(line,30)'IWETIT',this%Iwetit
      !  write(iu,60) trim(line)
      !endif
      !if (this%Ihdwet/=0) write(iu,30)'IHDWET', this%Ihdwet
    endif
    !!
    !if (.not. this%Newton) then
    !  write(iu,50)'NO_NEWTON'
    !endif
    !
!    if (this%BottomHeadDamping) then
!      write(iu,50)'BOTTOM_HEAD_DAMPENING'
!    endif
    !
    ! Write END Options
    write(iu,100)
    !
    return
  end subroutine WriteOptions

  double precision function get_q3d_hk() result(hk)
    ! Based on ITMUNI, LENUNI, and Domenico & Schwartz, return a
    ! reasonable, small value for horizontal hydraulic conductivity,
    ! to be used in active cells that will replace a quasi-3D
    ! confining unit.
    implicit none
    ! local
    ! -- hkms is a typical, small horizontal
    !    hydraulic conductivity in meters and seconds.
    double precision :: converter, hkms
    character(len=300) :: msg
    character(len=10)  :: strng
    ! format
    10 format(g10.3)
    !
    ! -- Start with a small HK value characteristic of shale
    !    or unfractured igneous or metamorphic rock
    !    (Domenico and Schwartz, 1990)
    hkms = 1.0e-13
    !
    ! -- Convert to model units specified in BAS input
    if (ITMUNI < 1 .or. ITMUNI > 5) then
      msg = 'Error: ITMUNI is undefined in Discretization Package input. ' // &
          'A valid value needs to be assigned so that HK can be ' // &
          'assigned for layer(s) representing former quasi-3D unit(s).'
      call store_error(msg)
    endif
    if (LENUNI < 1 .or. LENUNI > 3) then
      msg = 'Error: LENUNI is undefined in Discretization Package input. ' // &
          'A valid value needs to be assigned so that HK can be ' // &
          'assigned for layer(s) representing former quasi-3D unit(s).'
      call store_error(msg)
    endif
    if (count_errors() > 0) call ustop()
    !
    ! -- Define time conversion
    select case (ITMUNI)
    case (1)
      ! seconds -- no time conversion needed
      converter = DONE
    case (2)
      ! minutes
      converter = 60.0d0
    case (3)
      ! hours
      converter = 60.0d0 * 60.0d0
    case (4)
      ! days
      converter = 60.0d0 * 60.0d0 * 24.0d0
    case (5)
      ! years
      converter = 60.0d0 * 60.0d0 * 24.0d0 * 365.25d0
    end select
    !
    ! -- Define length conversion
    select case (LENUNI)
    case (1)
      ! feet
      converter = converter / 0.3048d0
    case (2)
      ! meters -- no length conversion needed
    case (3)
      ! centimeters
      converter = converter * 100.0d0
    end select
    !
    ! -- Perform conversion
    hk = hkms * converter
    ! -- Round the value to about 3 significant figures
    write(strng,10)hk
    read(strng,*)hk
    !
    return
  end function get_q3d_hk

  !
  ! Functions to assign VK based on Vcont and cell thicknesses from BCF
  !

  subroutine assign_vk(this, vkerr)
    ! Assign VK array for all layers when BCF is used and NLAY > 1,
    ! based on vcont arrays and cell thicknesses.
    implicit none
    ! dummy
    class(NpfWriterType) :: this
    logical, intent(inout) :: vkerr
    ! local
    integer :: i, j, k, knew
    double precision :: thk, thkp1, thq, vktemp
    double precision :: topk, botk, botq, botkp1
    double precision, parameter :: vkbig = 1.0d20
    character(len=1000) :: msg
    !
    if (.not. associated(vcont)) return
    vkerr = .false.
    knew = 0
    bigloop: do k=1,nlay-1
      knew = knew + 1
      if (LAYCBD(k) == 0) then
        ! No quasi-3d unit underlies layer k
        if (k == 1) then
          ! Assign VK for layers 1 and 2
          do i=1,NROW
            do j=1,NCOL
              if (k == 1) then
                topk = BOTM(j,i,0)
              else
                topk = BOTM(j,i,LBOTM(k-1))
              endif
              botk = BOTM(j,i,LBOTM(k))
              botkp1 = BOTM(j,i,LBOTM(k)+1)
              thk = topk - botk
              thkp1 = botk - botkp1
              vktemp = vk_noq3d(vcont(j,i,k), thk, thkp1)
              this%vk(j,i,knew) = vktemp
              this%vk(j,i,knew+1) = vktemp
            enddo
          enddo
        else
          ! VK for current layer has already been assigned.
          ! Assign VK for next layer.
          do i=1,NROW
            do j=1,NCOL
              topk = BOTM(j,i,LBOTM(k-1))
              botk = BOTM(j,i,LBOTM(k))
              botkp1 = BOTM(j,i,LBOTM(k)+1)
              thk = topk - botk
              thkp1 = botk - botkp1
              this%vk(j,i,knew+1) = vkkp1(vcont(j,i,k), thk, thkp1, &
                                        this%vk(j,i,knew), vkerr)
              if (vkerr) exit bigloop
            enddo
          enddo
        endif
      else
        ! A quasi-3d unit underlies layer k
        if (k == 1) then
          ! Assign VK for layers 1, 2, and 3
          do i=1,NROW
            do j=1,NCOL
              topk = BOTM(j,i,LBOTM(k)-1)
              botk = BOTM(j,i,LBOTM(k))
              botq = BOTM(j,i,LBOTM(k)+1)
              botkp1 = BOTM(j,i,LBOTM(k+1))
              thk = topk - botk
              thq = botk - botq
              thkp1 = botq - botkp1
              this%vk(j,i,knew) = vkbig
              this%vk(j,i,knew+1) = vcont(j,i,k) * thq
              this%vk(j,i,knew+2) = vkbig
            enddo
          enddo
        else
          ! VK for current layer has already been assigned.
          ! Assign VK for quasi-3d unit and next layer.
          do i=1,NROW
            do j=1,NCOL
              topk = BOTM(j,i,LBOTM(k)-1)
              botk = BOTM(j,i,LBOTM(k))
              botq = BOTM(j,i,LBOTM(k)+1)
              botkp1 = BOTM(j,i,LBOTM(k+1))
              thk = topk - botk
              thq = botk - botq
              thkp1 = botq - botkp1
              vktemp = vkcj(vcont(j,i,k), thk, thkp1, thq, &
                            this%vk(j,i,knew), vkerr)
              this%vk(j,i,knew+1) = vktemp
              this%vk(j,i,knew+2) = vktemp
            enddo
          enddo
        endif
        ! Increment knew
        knew = knew + 1
      endif
    enddo bigloop
    if (vkerr) then
      msg = 'Unable to assign physically meaningful VK values that' &
            // ' produce specified Vcont values.'
      call store_warning(msg)
      ! Ned todo: write and call a routine to assign placeholder values for
      ! VK arrays.
    else
      msg = 'VK values that are equivalent to specified Vcont values' &
            // ' have been assigned.  However, these values may or' &
            // ' may not be reasonable.'
      call store_note(msg)
    endif
    !
    return
  end subroutine assign_vk

  double precision function vk_noq3d(vcnt, thk, thkp1)
    ! For case without quasi-3d unit, vkk and vkkp1 unknown.
    ! Assume vkk = vkkp1 = vk.  Solve for vk.
    implicit none
    ! dummy
    double precision, intent(in) :: vcnt, thk, thkp1
    !
    vk_noq3d = vcnt*DHALF*(thk + thkp1)
    return
  end function vk_noq3d

  double precision function vkkp1(vcnt, thk, thkp1, vkk, vkerr)
    ! For case without quasi-3d unit, vkk known.
    ! Solve for vkkp1.
    implicit none
    ! dummy
    double precision, intent(in) :: vcnt, thk, thkp1, vkk
    logical, intent(inout) :: vkerr
    ! local
    double precision :: denom
    !
    denom = (DONE/vcnt - DHALF*thk/vkk)
    if (denom <= DZERO) then
      vkerr = .true.
      vkkp1 = -9999.0d0
    else
      vkkp1 = DHALF*thkp1/denom
    endif
    return
  end function vkkp1

  double precision function vkcj(vcnt, thk, thkp1, thq, vkk, vkerr)
    ! For case with quasi-3d unit, vkk known.
    ! Assume VKCB = vkkp1 = vkcj.  Solve for vkcj.
    implicit none
    ! dummy
    double precision, intent(in) :: vcnt, thk, thkp1, thq, vkk
    logical, intent(inout) :: vkerr
    ! local
    double precision :: denom
    !
    denom = DONE/vcnt - DHALF*thk/vkk
    if (denom <= DZERO) then
      vkerr = .true.
      vkcj = -9999.0d0
    else
      vkcj = (thq + DHALF*thkp1) / denom
    endif
    return
  end function vkcj

end module NpfWriterModule
