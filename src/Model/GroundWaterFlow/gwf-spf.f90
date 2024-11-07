module SpfModule
  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DZERO, DHALF, DEM6, LENFTYPE, LENPACKAGENAME
  use SimVariablesModule, only: errmsg
  use SimModule, only: count_errors, store_error, store_error_filename
  use MemoryManagerModule, only: mem_allocate, mem_deallocate, &
                                 mem_setptr, mem_checkin
  use MemoryHelperModule, only: create_mem_path
  use BndModule, only: BndType
  use BndExtModule, only: BndExtType
  use MatrixBaseModule
  
  implicit none
  private

  public :: spf_create
  public :: SpfType
  
  character(len=LENFTYPE) :: ftype = 'SPF'
  character(len=LENPACKAGENAME) :: text = '             SPF'
  !
  type, extends(BndExtType) :: SpfType
    real(DP), dimension(:), pointer, contiguous :: distance => null() !< SPF distance to face
    real(DP), dimension(:), pointer, contiguous :: area => null() !< SPF face area
    logical(LGP), private, pointer :: some_option => null() !< some option
  contains
    procedure :: allocate_scalars => spf_allocate_scalars
    procedure :: allocate_arrays => spf_allocate_arrays
    procedure :: source_options => spf_source_options
    procedure :: bnd_cf => spf_cf
    procedure :: bnd_fc => spf_fc
    procedure :: bnd_da => spf_da
    procedure :: define_listlabel
  end type SpfType

contains

  !> @brief Create a New SPF Package
  !<
  subroutine spf_create(packobj, id, ibcnum, inunit, iout, namemodel, pakname, &
                        mempath)
    class(BndType), pointer :: packobj
    integer(I4B), intent(in) :: id
    integer(I4B), intent(in) :: ibcnum
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    character(len=*), intent(in) :: namemodel
    character(len=*), intent(in) :: pakname
    character(len=*), intent(in) :: mempath
    ! local
    type(SpfType), pointer :: spfobj
    
    ! allocate the object and assign values to object variables
    allocate (spfobj)
    packobj => spfobj
    
    ! create name and memory path
    call packobj%set_names(ibcnum, namemodel, pakname, ftype, mempath)
    packobj%text = text
    
    ! allocate scalars
    call packobj%allocate_scalars()
    
    ! initialize package
    call packobj%pack_initialize()
    
    packobj%inunit = inunit
    packobj%iout = iout
    packobj%id = id
    packobj%ibcnum = ibcnum
    packobj%ictMemPath = create_mem_path(namemodel, 'NPF') ! TOD_UZR: why do we have this?

    spfobj%some_option = .false.
    
  end subroutine spf_create

  subroutine spf_source_options(this)
    use MemoryManagerExtModule, only: mem_set_value
    use GwfSpfInputModule, only: GwfSpfParamFoundType
    class(SpfType), intent(inout) :: this
    ! local
    type(GwfSpfParamFoundType) :: found
    
    ! source common bound options
    call this%BndExtType%source_options()

    call mem_set_value(this%some_option, 'SOME_OPTION', &
                       this%input_mempath, found%some_option)

  end subroutine spf_source_options

  !> @brief Allocate scalars
  !<
  subroutine spf_allocate_scalars(this)
    class(SpfType) :: this !< this instance

    ! base allocate
    call this%BndExtType%allocate_scalars()

    call mem_allocate(this%some_option, 'DEV_VS2D_BND', this%memoryPath)

  end subroutine spf_allocate_scalars

  !> @brief Allocate arrays
  !<
  subroutine spf_allocate_arrays(this, nodelist, auxvar)
    class(SpfType) :: this
    integer(I4B), dimension(:), pointer, contiguous, optional :: nodelist
    real(DP), dimension(:, :), pointer, contiguous, optional :: auxvar
    
    ! call base type allocate arrays
    call this%BndExtType%allocate_arrays(nodelist, auxvar)
    
    ! set spf input context pointers
    call mem_setptr(this%distance, 'DIST', this%input_mempath)
    call mem_setptr(this%area, 'AREA', this%input_mempath)
    
    ! checkin spf input context pointers
    call mem_checkin(this%distance, 'FACEDIST', this%memoryPath, &
                     'DIST', this%input_mempath)
    call mem_checkin(this%area, 'FACEAREA', this%memoryPath, &
                     'AREA', this%input_mempath)
  
  end subroutine spf_allocate_arrays

  !> @brief Formulate the HCOF and RHS terms
  !<
  subroutine spf_cf(this)
    class(SpfType) :: this
    ! local
    integer(I4B) :: i, node
    real(DP) :: z, head
    real(DP), dimension(:), contiguous, pointer :: condsat
    
    if (this%nbound .eq. 0) return

    call mem_setptr(condsat, 'CONDSAT', this%ictMemPath) ! TODO_UZR: remove this HACK
    
    ! Calculate hcof and rhs for each seepage face
    do i = 1, this%nbound
      node = this%nodelist(i)
      if (this%ibound(node) <= 0) then
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
        cycle
      end if

      z = DHALF * (this%dis%bot(node) + this%dis%top(node))
      head = this%xnew(node)
      if (head > z) then
        this%hcof(i) = -this%area(i) * condsat(node) / this%distance(i)
        this%rhs(i) = -this%area(i) * condsat(node) * z / this%distance(i)
      else
        this%hcof(i) = DZERO
        this%rhs(i) = DZERO
      end if
      
    end do
    
  end subroutine spf_cf

  !> @brief Copy rhs and hcof into solution rhs and amat
  !<
  subroutine spf_fc(this, rhs, ia, idxglo, matrix_sln)
    class(SpfType) :: this
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B), dimension(:), intent(in) :: ia
    integer(I4B), dimension(:), intent(in) :: idxglo
    class(MatrixBaseType), pointer :: matrix_sln
    ! local
    integer(I4B) :: i, n, ipos

    ! pakmvrobj fc
    if (this%imover == 1) then
      call this%pakmvrobj%fc()
    end if
    
    ! Copy package rhs and hcof into solution rhs and amat
    do i = 1, this%nbound
      n = this%nodelist(i)
      rhs(n) = rhs(n) + this%rhs(i)
      ipos = ia(n)
      call matrix_sln%add_value_pos(idxglo(ipos), this%hcof(i))
      
      ! If mover is active and this boundary is discharging,
      ! store available water (as positive value).
      ! TODO_UZR: implement mover
      if (this%imover == 1) then
        write(errmsg, '(a,a)') "Seepage Mover not supported: ", this%packName
        call store_error(errmsg, terminate=.true.)
      end if
    end do
    
  end subroutine spf_fc

  !> @brief Define the list heading that is written to iout when PRINT_INPUT
  !< option is used
  subroutine define_listlabel(this)
    
    class(SpfType), intent(inout) :: this
    
    ! create the header list label
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
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'FACE DISTANCE'
    write (this%listlabel, '(a, a16)') trim(this%listlabel), 'FACE AREA'
    if (this%inamedbound == 1) then
      write (this%listlabel, '(a, a16)') trim(this%listlabel), 'BOUNDARY NAME'
    end if
    
  end subroutine define_listlabel
  
  !> @brief Deallocate memory
  !<
  subroutine spf_da(this)
    class(SpfType) :: this
    
    call this%BndExtType%bnd_da()
    
    call mem_deallocate(this%distance, 'FACEDIST', this%memoryPath)
    call mem_deallocate(this%area, 'FACEAREA', this%memoryPath)

    call mem_deallocate(this%some_option)
   
  end subroutine spf_da

end module SpfModule
